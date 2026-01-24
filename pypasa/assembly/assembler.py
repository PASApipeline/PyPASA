"""
PASA assembly wrapper for PyPASA.

This module wraps the existing C++ PASA assembly utility (pasa_cpp/pasa)
which contains the proven assembly algorithm. Rather than rewriting the
complex assembly logic in Python, we:
1. Export clusters to temp files
2. Call the C++ binary
3. Parse the assembly output
4. Load assemblies into the DataStore

Corresponds to assemble_clusters.dbi and assembly_db_loader.dbi in Perl PASA.
"""

import os
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import List, Dict, Optional
import pickle

from pypasa.storage.models import Alignment, Assembly, AlignmentSegment, Orientation
from pypasa.storage.datastore import DataStore

logger = logging.getLogger(__name__)


class PASAAssembler:
    """
    Wrapper for C++ PASA assembly utility.
    
    The assembly algorithm remains in C++ (pasa binary from pasa_cpp/).
    This class handles data export, execution, and result import.
    """
    
    def __init__(
        self,
        datastore: DataStore,
        pasa_home: str,
        genome_fasta: str,
        fuzzlength: int = 20,
    ):
        """
        Initialize assembler.
        
        Args:
            datastore: DataStore for querying/storing data
            pasa_home: Path to PASA installation (contains bin/pasa)
            genome_fasta: Path to genome FASTA file
            fuzzlength: Fuzzy length for assembly (default 20)
        """
        self.datastore = datastore
        self.pasa_home = Path(pasa_home)
        self.genome_fasta = Path(genome_fasta)
        self.fuzzlength = fuzzlength
        
        # Locate C++ pasa binary
        self.pasa_bin = self.pasa_home / "bin" / "pasa"
        if not self.pasa_bin.exists():
            raise FileNotFoundError(
                f"PASA binary not found: {self.pasa_bin}. "
                f"Have you run 'make' in {self.pasa_home}?"
            )
        
        self.next_assembly_id = 1
        self.stats = {
            'clusters_processed': 0,
            'assemblies_generated': 0,
            'failed_clusters': 0,
        }
    
    def assemble_all_clusters(
        self, 
        output_dir: str = "./assemblies",
        scaffold: str = None
    ) -> Dict[str, int]:
        """
        Assemble all clusters.
        
        Args:
            output_dir: Directory to store assembly files
            scaffold: If provided, only assemble clusters on this scaffold
        
        Returns:
            Dictionary of assembly statistics
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Starting assembly of all clusters")
        
        # Get scaffolds to process
        if scaffold:
            scaffolds = [scaffold]
        else:
            scaffolds = self.datastore.get_scaffolds()
        
        for scaffold in scaffolds:
            self._assemble_scaffold_clusters(scaffold, output_dir)
        
        logger.info(
            f"Assembly complete: {self.stats['assemblies_generated']} assemblies "
            f"from {self.stats['clusters_processed']} clusters"
        )
        
        return self.stats
    
    def _assemble_scaffold_clusters(self, scaffold: str, output_dir: Path) -> None:
        """
        Assemble all clusters on a scaffold.
        
        Args:
            scaffold: Scaffold name
            output_dir: Output directory
        """
        # Load alignments for this scaffold
        alignments = self.datastore.get_alignments_by_scaffold(scaffold, valid_only=True)
        
        if not alignments:
            return
        
        # Group alignments by cluster
        clusters = self._group_by_cluster(alignments)
        
        logger.debug(f"Assembling {len(clusters)} clusters on {scaffold}")
        
        for cluster_id, cluster_alignments in clusters.items():
            try:
                assemblies = self._assemble_cluster(
                    cluster_id, cluster_alignments, scaffold, output_dir
                )
                
                if assemblies:
                    self.datastore.add_assemblies(assemblies)
                    self.stats['assemblies_generated'] += len(assemblies)
                
                self.stats['clusters_processed'] += 1
                
            except Exception as e:
                logger.warning(f"Failed to assemble cluster {cluster_id}: {e}")
                self.stats['failed_clusters'] += 1
                continue
    
    def _group_by_cluster(
        self, alignments: List[Alignment]
    ) -> Dict[int, List[Alignment]]:
        """
        Group alignments by cluster ID.
        
        Args:
            alignments: List of alignments
            
        Returns:
            Dictionary mapping cluster_id -> List[Alignment]
        """
        from collections import defaultdict
        
        clusters = defaultdict(list)
        for alignment in alignments:
            if alignment.cluster_id is not None:
                clusters[alignment.cluster_id].append(alignment)
        
        return dict(clusters)
    
    def _assemble_cluster(
        self,
        cluster_id: int,
        alignments: List[Alignment],
        scaffold: str,
        output_dir: Path,
    ) -> List[Assembly]:
        """
        Assemble a single cluster using C++ pasa utility.
        
        Args:
            cluster_id: Cluster ID
            alignments: Alignments in this cluster
            scaffold: Scaffold name
            output_dir: Output directory
            
        Returns:
            List of Assembly objects generated
        """
        if not alignments:
            return []
        
        # Create temp file for cluster data
        with tempfile.NamedTemporaryFile(
            mode='wb', suffix='.cluster', delete=False, dir=output_dir
        ) as cluster_file:
            cluster_path = cluster_file.name
            # Serialize alignments (Python pickle - the C++ tool will need adjustments)
            # For now, we'll use a text format that matches what assemble_clusters.dbi exports
            pickle.dump(alignments, cluster_file)
        
        try:
            # For now, create a placeholder assembly since integrating with C++ requires
            # matching the data format it expects. We'll implement proper integration later.
            assemblies = self._create_placeholder_assembly(cluster_id, alignments, scaffold)
            
            return assemblies
            
        finally:
            # Clean up temp file
            if Path(cluster_path).exists():
                Path(cluster_path).unlink()
    
    def _create_placeholder_assembly(
        self,
        cluster_id: int,
        alignments: List[Alignment],
        scaffold: str,
    ) -> List[Assembly]:
        """
        Create a simple assembly from alignments.
        
        This is a placeholder until we properly integrate with the C++ pasa utility.
        For now, we just take the longest alignment as the assembly.
        
        Args:
            cluster_id: Cluster ID
            alignments: Alignments to assemble
            scaffold: Scaffold name
            
        Returns:
            List containing single Assembly
        """
        # Sort by length and take longest
        longest = max(alignments, key=lambda a: a.aligned_length)
        
        # Create assembly from longest alignment
        assembly = Assembly(
            assembly_id=self.next_assembly_id,
            cluster_id=cluster_id,
            scaffold=scaffold,
            lend=longest.genome_lend,
            rend=longest.genome_rend,
            orient=longest.orient,
            segments=longest.segments.copy(),
            component_align_ids=[a.align_id for a in alignments],
            sequence=None,  # Would need to extract from genome
            is_fli=longest.is_fli,
        )
        
        self.next_assembly_id += 1
        
        return [assembly]
    
    def _call_pasa_binary(
        self,
        cluster_file: str,
        scaffold: str,
        output_file: str,
    ) -> subprocess.CompletedProcess:
        """
        Call the C++ pasa binary.
        
        Args:
            cluster_file: Path to serialized cluster file
            scaffold: Scaffold name
            output_file: Path to output file
            
        Returns:
            CompletedProcess result
        """
        cmd = [
            str(self.pasa_bin),
            "-f", str(self.fuzzlength),
            "-i", cluster_file,
            "-s", scaffold,
            "-g", str(self.genome_fasta),
            "-o", output_file,
        ]
        
        logger.debug(f"Running: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        
        return result
    
    def _parse_assembly_output(self, output_file: str) -> List[Assembly]:
        """
        Parse assembly output from C++ pasa utility.
        
        Args:
            output_file: Path to assembly output file
            
        Returns:
            List of Assembly objects
        """
        # TODO: Implement proper parsing of pasa output format
        # For now, return empty list
        return []
    
    def get_assembly_report(self) -> str:
        """
        Generate a human-readable assembly report.
        
        Returns:
            Formatted report string
        """
        report = []
        report.append("=" * 60)
        report.append("PASA ASSEMBLY REPORT")
        report.append("=" * 60)
        report.append(f"Clusters processed:      {self.stats['clusters_processed']:>10,}")
        report.append(f"Assemblies generated:    {self.stats['assemblies_generated']:>10,}")
        report.append(f"Failed clusters:         {self.stats['failed_clusters']:>10,}")
        
        if self.stats['clusters_processed'] > 0:
            avg = self.stats['assemblies_generated'] / self.stats['clusters_processed']
            report.append(f"Avg assemblies/cluster:  {avg:>10.2f}")
        
        report.append("=" * 60)
        
        return "\n".join(report)
