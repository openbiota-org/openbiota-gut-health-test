"""Expanded gut bacterial and archaeal recognition (spec 0.8.4).

Modules
-------
gtdb        GTDB R232 backbone: slim genome table, species clusters, radii.
crosswalk   Release-aware SGB -> GTDB crosswalks built from assembly
            membership, with typed relationships (exact, synonym, split,
            merge, unresolved). Never a numeric-suffix join.
globdb      GlobDB r232 taxonomy and source dictionaries.
reconcile   R5-R10 supplements against GTDB232/GlobDB: already_represented,
            new_cluster, additional_strain_reference, better_reference,
            unresolved_mapping, excluded_quality; the rescue-panel manifest.
locks       ReferenceRelease records from refs/expanded/locks and the tool lock.

Everything here produces files under refs/ once and records how; sample
runs read those files and never resolve a moving alias.
"""
