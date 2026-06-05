=============
API reference
=============

This page documents the public Cell-GPS functions intended for direct use in
analysis scripts. The preferred import path is ``cellgps``; the legacy
``sfplot`` namespace remains available for backward compatibility.

Core COSTE and StructureMap functions
-------------------------------------

.. autofunction:: cellgps.compute_cophenetic_distances_from_df

.. autofunction:: cellgps.compute_cophenetic_distances_from_adata

.. autofunction:: cellgps.compute_searcher_findee_distance_matrix_from_df

.. autofunction:: cellgps.compute_cophenetic_from_distance_matrix

.. autofunction:: cellgps.compute_cophenetic_distances_from_df_memory_opt

.. autofunction:: cellgps.pick_batch_size

Topology extensions
-------------------

.. autofunction:: cellgps.compute_weighted_cophenetic_distances_from_df

.. autofunction:: cellgps.compute_weighted_searcher_findee_distance_matrix_from_df

.. autofunction:: cellgps.build_entity_points_from_expression

.. autofunction:: cellgps.compute_entity_structuremap

.. autofunction:: cellgps.compute_entity_to_cell_topology

.. autofunction:: cellgps.compute_pathway_activity_matrix

.. autofunction:: cellgps.ligand_receptor_topology_analysis

.. autofunction:: cellgps.ligand_receptor_target_consistency

.. autofunction:: cellgps.pathway_topology_analysis

Preprocessing and input helpers
-------------------------------

.. autofunction:: cellgps.load_xenium_data

.. autofunction:: cellgps.load_xenium_table_bundle

.. autofunction:: cellgps.merge_xenium_clusters_into_adata

.. autofunction:: cellgps.read_visium_bin

Plotting
--------

.. autofunction:: cellgps.plot_cophenetic_heatmap

.. autofunction:: cellgps.generate_cluster_distance_heatmap_from_adata

.. autofunction:: cellgps.generate_cluster_distance_heatmap_from_df

.. autofunction:: cellgps.generate_cluster_distance_heatmap_from_path

.. autofunction:: cellgps.plot_circular_dendrogram_pycirclize

.. autofunction:: cellgps.circle_heatmap
