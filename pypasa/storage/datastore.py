"""
Hybrid DataStore for PyPASA.

This module provides a storage manager that combines:
- SQLite for persistence and queryability
- In-memory structures for fast processing
- Checkpointing for recovery

Design:
- SQLite stores all data for durability
- Data loaded into RAM for processing
- Process one scaffold at a time to manage memory
- Write-ahead log (WAL) mode for concurrent access
"""

import sqlite3
import json
import pickle
from pathlib import Path
from typing import List, Optional, Dict, Iterator
from contextlib import contextmanager
import logging

from pypasa.storage.models import Alignment, Assembly, AlignmentSegment, Orientation

logger = logging.getLogger(__name__)


class DataStore:
    """
    Hybrid storage manager for PyPASA alignments and assemblies.
    
    Provides:
    - Persistent SQLite storage
    - Bulk load/save operations
    - Scaffold-level data access
    - Checkpoint/resume support
    """
    
    def __init__(self, db_path: str):
        """
        Initialize datastore.
        
        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.conn: Optional[sqlite3.Connection] = None
        self._init_database()
        
    def _init_database(self) -> None:
        """Initialize SQLite database with schema."""
        self.conn = sqlite3.connect(
            str(self.db_path),
            check_same_thread=False,  # Allow use across threads
            isolation_level=None,  # Autocommit mode
        )
        
        # Enable Write-Ahead Logging for better concurrency
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute("PRAGMA cache_size=-64000")  # 64MB cache
        
        # Create schema
        self.conn.executescript("""
            -- Alignments
            CREATE TABLE IF NOT EXISTS alignments (
                align_id INTEGER PRIMARY KEY,
                cdna_acc TEXT NOT NULL,
                scaffold TEXT NOT NULL,
                lend INTEGER NOT NULL,
                rend INTEGER NOT NULL,
                orient TEXT NOT NULL,
                avg_percent_id REAL,
                percent_aligned REAL,
                valid INTEGER DEFAULT 0,
                cluster_id INTEGER,
                is_fli INTEGER DEFAULT 0,
                spliced_orient TEXT,
                num_segments INTEGER,
                alignment_data TEXT  -- JSON serialized
            );
            
            CREATE INDEX IF NOT EXISTS idx_align_scaffold 
                ON alignments(scaffold, lend, rend);
            CREATE INDEX IF NOT EXISTS idx_align_cluster 
                ON alignments(cluster_id) WHERE cluster_id IS NOT NULL;
            CREATE INDEX IF NOT EXISTS idx_align_cdna 
                ON alignments(cdna_acc);
            CREATE INDEX IF NOT EXISTS idx_align_valid 
                ON alignments(valid) WHERE valid = 1;
            
            -- Alignment segments (for detailed structure)
            CREATE TABLE IF NOT EXISTS alignment_segments (
                align_id INTEGER NOT NULL,
                segment_idx INTEGER NOT NULL,
                genome_lend INTEGER NOT NULL,
                genome_rend INTEGER NOT NULL,
                transcript_lend INTEGER NOT NULL,
                transcript_rend INTEGER NOT NULL,
                percent_id REAL,
                PRIMARY KEY (align_id, segment_idx),
                FOREIGN KEY (align_id) REFERENCES alignments(align_id)
            );
            
            -- Assemblies
            CREATE TABLE IF NOT EXISTS assemblies (
                assembly_id INTEGER PRIMARY KEY,
                cluster_id INTEGER NOT NULL,
                scaffold TEXT NOT NULL,
                lend INTEGER NOT NULL,
                rend INTEGER NOT NULL,
                orient TEXT NOT NULL,
                num_segments INTEGER,
                is_fli INTEGER DEFAULT 0,
                sequence TEXT,
                assembly_data TEXT  -- JSON serialized
            );
            
            CREATE INDEX IF NOT EXISTS idx_assembly_scaffold 
                ON assemblies(scaffold, lend, rend);
            CREATE INDEX IF NOT EXISTS idx_assembly_cluster 
                ON assemblies(cluster_id);
            
            -- Assembly components (which alignments went into each assembly)
            CREATE TABLE IF NOT EXISTS assembly_components (
                assembly_id INTEGER NOT NULL,
                align_id INTEGER NOT NULL,
                PRIMARY KEY (assembly_id, align_id),
                FOREIGN KEY (assembly_id) REFERENCES assemblies(assembly_id),
                FOREIGN KEY (align_id) REFERENCES alignments(align_id)
            );
            
            -- Metadata for checkpointing
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT
            );
        """)
        
        logger.info(f"Initialized database: {self.db_path}")
    
    def close(self) -> None:
        """Close database connection."""
        if self.conn:
            self.conn.close()
            logger.info("Closed database connection")
    
    @contextmanager
    def transaction(self):
        """Context manager for explicit transactions."""
        cursor = self.conn.cursor()
        try:
            cursor.execute("BEGIN")
            yield cursor
            cursor.execute("COMMIT")
        except Exception:
            cursor.execute("ROLLBACK")
            raise
    
    # ========== Alignment Operations ==========
    
    def add_alignments(self, alignments: List[Alignment]) -> None:
        """
        Bulk insert alignments.
        
        Args:
            alignments: List of Alignment objects
        """
        if not alignments:
            return
        
        with self.transaction() as cursor:
            # Insert alignment records
            cursor.executemany(
                """
                INSERT INTO alignments (
                    align_id, cdna_acc, scaffold, lend, rend, orient,
                    avg_percent_id, percent_aligned, valid, cluster_id,
                    is_fli, spliced_orient, num_segments, alignment_data
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        aln.align_id,
                        aln.cdna_acc,
                        aln.scaffold,
                        aln.genome_lend,
                        aln.genome_rend,
                        aln.orient.value,
                        aln.avg_percent_id,
                        aln.percent_aligned,
                        1 if aln.valid else 0,
                        aln.cluster_id,
                        1 if aln.is_fli else 0,
                        aln.spliced_orient.value if aln.spliced_orient else None,
                        aln.num_segments,
                        json.dumps(aln.to_dict()),
                    )
                    for aln in alignments
                ],
            )
            
            # Insert segment records
            segment_rows = []
            for aln in alignments:
                for idx, seg in enumerate(aln.segments):
                    segment_rows.append(
                        (
                            aln.align_id,
                            idx,
                            seg.genome_lend,
                            seg.genome_rend,
                            seg.transcript_lend,
                            seg.transcript_rend,
                            seg.percent_id,
                        )
                    )
            
            if segment_rows:
                cursor.executemany(
                    """
                    INSERT INTO alignment_segments (
                        align_id, segment_idx, genome_lend, genome_rend,
                        transcript_lend, transcript_rend, percent_id
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    segment_rows,
                )
        
        logger.info(f"Added {len(alignments)} alignments")
    
    def get_alignments_by_scaffold(
        self, scaffold: str, valid_only: bool = False
    ) -> List[Alignment]:
        """
        Load all alignments for a scaffold.
        
        Args:
            scaffold: Scaffold name
            valid_only: If True, only return valid alignments
            
        Returns:
            List of Alignment objects
        """
        query = "SELECT alignment_data FROM alignments WHERE scaffold = ?"
        params = [scaffold]
        
        if valid_only:
            query += " AND valid = 1"
        
        cursor = self.conn.execute(query, params)
        
        alignments = []
        for (data,) in cursor:
            aln_dict = json.loads(data)
            alignments.append(self._dict_to_alignment(aln_dict))
        
        logger.debug(f"Loaded {len(alignments)} alignments for scaffold {scaffold}")
        return alignments
    
    def update_alignment_clusters(self, cluster_assignments: Dict[int, int]) -> None:
        """
        Update cluster_id for alignments.
        
        Args:
            cluster_assignments: Mapping of align_id -> cluster_id
        """
        if not cluster_assignments:
            return
        
        with self.transaction() as cursor:
            cursor.executemany(
                "UPDATE alignments SET cluster_id = ? WHERE align_id = ?",
                [(cluster_id, align_id) for align_id, cluster_id in cluster_assignments.items()],
            )
        
        logger.info(f"Updated cluster assignments for {len(cluster_assignments)} alignments")
    
    def update_alignment_validity(self, validity: Dict[int, bool]) -> None:
        """
        Update valid flag for alignments.
        
        Args:
            validity: Mapping of align_id -> valid (True/False)
        """
        if not validity:
            return
        
        with self.transaction() as cursor:
            # Update the valid column
            cursor.executemany(
                "UPDATE alignments SET valid = ? WHERE align_id = ?",
                [(1 if is_valid else 0, align_id) for align_id, is_valid in validity.items()],
            )
            
            # Also update the JSON data to keep it in sync
            for align_id, is_valid in validity.items():
                cursor.execute(
                    """
                    UPDATE alignments 
                    SET alignment_data = json_set(alignment_data, '$.valid', ?)
                    WHERE align_id = ?
                    """,
                    (1 if is_valid else 0, align_id),
                )
        
        logger.info(f"Updated validity for {len(validity)} alignments")
    
    def update_percent_aligned(self, percent_aligned: Dict[int, float]) -> None:
        """
        Update percent_aligned values for alignments.
        
        Args:
            percent_aligned: Mapping of align_id -> percent_aligned
        """
        if not percent_aligned:
            return
        
        with self.transaction() as cursor:
            # Update the percent_aligned column
            cursor.executemany(
                "UPDATE alignments SET percent_aligned = ? WHERE align_id = ?",
                [(pct, align_id) for align_id, pct in percent_aligned.items()],
            )
            
            # Also update the JSON data to keep it in sync
            for align_id, pct in percent_aligned.items():
                cursor.execute(
                    """
                    UPDATE alignments 
                    SET alignment_data = json_set(alignment_data, '$.percent_aligned', ?)
                    WHERE align_id = ?
                    """,
                    (pct, align_id),
                )
        
        logger.debug(f"Updated percent_aligned for {len(percent_aligned)} alignments")
    
    def update_spliced_orientations(self, orientations: Dict[int, 'Orientation']) -> None:
        """
        Update spliced orientation (transcribed orientation) for alignments.
        
        Args:
            orientations: Mapping of align_id -> Orientation
        """
        if not orientations:
            return
        
        with self.transaction() as cursor:
            # Update the spliced_orient column
            cursor.executemany(
                "UPDATE alignments SET spliced_orient = ? WHERE align_id = ?",
                [(orient.value, align_id) for align_id, orient in orientations.items()],
            )
            
            # Also update the JSON data to keep it in sync
            for align_id, orient in orientations.items():
                cursor.execute(
                    """
                    UPDATE alignments 
                    SET alignment_data = json_set(alignment_data, '$.spliced_orient', ?)
                    WHERE align_id = ?
                    """,
                    (orient.value, align_id),
                )
        
        logger.info(f"Updated spliced orientation for {len(orientations)} alignments")
    
    def invalidate_alignments(self, align_ids: List[int]) -> None:
        """
        Mark alignments as invalid (used for best-alignment selection).
        
        Args:
            align_ids: List of alignment IDs to invalidate
        """
        if not align_ids:
            return
        
        with self.transaction() as cursor:
            cursor.executemany(
                "UPDATE alignments SET valid = 0 WHERE align_id = ?",
                [(align_id,) for align_id in align_ids],
            )
            
            # Also update the JSON data to keep it in sync
            cursor.executemany(
                """
                UPDATE alignments 
                SET alignment_data = json_set(alignment_data, '$.valid', 0)
                WHERE align_id = ?
                """,
                [(align_id,) for align_id in align_ids],
            )
        
        logger.info(f"Invalidated {len(align_ids)} alignments")
    
    def get_scaffolds(self) -> List[str]:
        """Get list of all scaffolds with alignments."""
        cursor = self.conn.execute(
            "SELECT DISTINCT scaffold FROM alignments ORDER BY scaffold"
        )
        return [row[0] for row in cursor]
    
    def get_alignments_by_cluster(self, cluster_id: int) -> List[Alignment]:
        """
        Get all alignments for a specific cluster.
        
        Args:
            cluster_id: Cluster ID to retrieve alignments for
            
        Returns:
            List of Alignment objects in this cluster
        """
        cursor = self.conn.execute(
            """
            SELECT align_id, alignment_data
            FROM alignments
            WHERE cluster_id = ?
            ORDER BY scaffold, lend
            """,
            (cluster_id,),
        )
        
        alignments = []
        for row in cursor:
            align_id, json_data = row
            data = json.loads(json_data)
            # Reconstruct Alignment from JSON
            segments = [AlignmentSegment(**seg) for seg in data['segments']]
            aln = Alignment(
                align_id=data['align_id'],
                cdna_acc=data['cdna_acc'],
                scaffold=data['scaffold'],
                orient=data['orient'],
                segments=segments,
                valid=bool(data.get('valid', 1)),
                avg_percent_id=data.get('avg_percent_id', 0.0),
                percent_aligned=data.get('percent_aligned', 0.0),
                cluster_id=data.get('cluster_id'),
            )
            alignments.append(aln)
        
        return alignments
    
    def get_clusters(self) -> List:
        """
        Get list of all clusters.
        
        Returns:
            List of namedtuples with cluster_id and alignment_count
        """
        from collections import namedtuple
        
        Cluster = namedtuple('Cluster', ['cluster_id', 'alignment_count'])
        
        cursor = self.conn.execute(
            """
            SELECT cluster_id, COUNT(*) as alignment_count
            FROM alignments
            WHERE cluster_id IS NOT NULL
            GROUP BY cluster_id
            ORDER BY cluster_id
            """
        )
        
        return [Cluster(cluster_id=row[0], alignment_count=row[1]) for row in cursor]
    
    # ========== Assembly Operations ==========
    
    def add_assemblies(self, assemblies: List[Assembly]) -> None:
        """
        Bulk insert assemblies.
        
        Args:
            assemblies: List of Assembly objects
        """
        if not assemblies:
            return
        
        with self.transaction() as cursor:
            # Insert assembly records
            cursor.executemany(
                """
                INSERT INTO assemblies (
                    assembly_id, cluster_id, scaffold, lend, rend, orient,
                    num_segments, is_fli, sequence, assembly_data
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        asm.assembly_id,
                        asm.cluster_id,
                        asm.scaffold,
                        asm.lend,
                        asm.rend,
                        asm.orient.value,
                        asm.num_segments,
                        1 if asm.is_fli else 0,
                        asm.sequence,
                        json.dumps(asm.to_dict()),
                    )
                    for asm in assemblies
                ],
            )
            
            # Insert component relationships
            component_rows = []
            for asm in assemblies:
                for align_id in asm.component_align_ids:
                    component_rows.append((asm.assembly_id, align_id))
            
            if component_rows:
                cursor.executemany(
                    "INSERT INTO assembly_components (assembly_id, align_id) VALUES (?, ?)",
                    component_rows,
                )
        
        logger.info(f"Added {len(assemblies)} assemblies")
    
    def get_assemblies_by_scaffold(self, scaffold: str) -> List[Assembly]:
        """Load all assemblies for a scaffold."""
        cursor = self.conn.execute(
            "SELECT assembly_data FROM assemblies WHERE scaffold = ?", (scaffold,)
        )
        
        assemblies = []
        for (data,) in cursor:
            asm_dict = json.loads(data)
            assemblies.append(self._dict_to_assembly(asm_dict))
        
        return assemblies
    
    # ========== Statistics ==========
    
    def get_stats(self) -> Dict[str, int]:
        """Get database statistics."""
        stats = {}
        
        cursor = self.conn.execute("SELECT COUNT(*) FROM alignments")
        stats["total_alignments"] = cursor.fetchone()[0]
        
        cursor = self.conn.execute("SELECT COUNT(*) FROM alignments WHERE valid = 1")
        stats["valid_alignments"] = cursor.fetchone()[0]
        
        cursor = self.conn.execute("SELECT COUNT(DISTINCT cluster_id) FROM alignments WHERE cluster_id IS NOT NULL")
        stats["clusters"] = cursor.fetchone()[0]
        
        cursor = self.conn.execute("SELECT COUNT(*) FROM assemblies")
        stats["assemblies"] = cursor.fetchone()[0]
        
        cursor = self.conn.execute("SELECT COUNT(DISTINCT scaffold) FROM alignments")
        stats["scaffolds"] = cursor.fetchone()[0]
        
        return stats
    
    # ========== Checkpointing ==========
    
    def set_metadata(self, key: str, value: str) -> None:
        """Store metadata key-value pair."""
        self.conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)", (key, value)
        )
    
    def get_metadata(self, key: str) -> Optional[str]:
        """Retrieve metadata value."""
        cursor = self.conn.execute("SELECT value FROM metadata WHERE key = ?", (key,))
        row = cursor.fetchone()
        return row[0] if row else None
    
    def checkpoint(self, stage: str) -> None:
        """Mark a processing stage as complete."""
        self.set_metadata(f"checkpoint_{stage}", "complete")
        logger.info(f"Checkpoint: {stage} complete")
    
    def is_checkpoint_complete(self, stage: str) -> bool:
        """Check if a checkpoint stage is complete."""
        return self.get_metadata(f"checkpoint_{stage}") == "complete"
    
    # ========== Helpers ==========
    
    @staticmethod
    def _dict_to_alignment(d: dict) -> Alignment:
        """Convert dictionary to Alignment object."""
        segments = [
            AlignmentSegment(
                genome_lend=s["genome_lend"],
                genome_rend=s["genome_rend"],
                transcript_lend=s["transcript_lend"],
                transcript_rend=s["transcript_rend"],
                percent_id=s["percent_id"],
            )
            for s in d["segments"]
        ]
        
        return Alignment(
            align_id=d["align_id"],
            cdna_acc=d["cdna_acc"],
            scaffold=d["scaffold"],
            segments=segments,
            orient=Orientation(d["orient"]),
            spliced_orient=Orientation(d["spliced_orient"]) if d.get("spliced_orient") else None,
            valid=d["valid"],
            cluster_id=d.get("cluster_id"),
            avg_percent_id=d.get("avg_percent_id"),
            percent_aligned=d.get("percent_aligned"),
            is_fli=d.get("is_fli", False),
        )
    
    @staticmethod
    def _dict_to_assembly(d: dict) -> Assembly:
        """Convert dictionary to Assembly object."""
        segments = [
            AlignmentSegment(
                genome_lend=s["genome_lend"],
                genome_rend=s["genome_rend"],
                transcript_lend=s["transcript_lend"],
                transcript_rend=s["transcript_rend"],
                percent_id=s["percent_id"],
            )
            for s in d["segments"]
        ]
        
        return Assembly(
            assembly_id=d["assembly_id"],
            cluster_id=d["cluster_id"],
            scaffold=d["scaffold"],
            lend=d["lend"],
            rend=d["rend"],
            orient=Orientation(d["orient"]),
            segments=segments,
            component_align_ids=d.get("component_align_ids", []),
            sequence=d.get("sequence"),
            is_fli=d.get("is_fli", False),
        )
