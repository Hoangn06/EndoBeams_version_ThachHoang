module EndoBeams

    include("packages.jl")
    include("constants.jl")
    include("linearsolvers.jl")
    include("def_params.jl") 

    include("def_bcs_and_loads.jl")
    include("def_configurations.jl")
    include("def_structures/def_structures_superelastic.jl")
    include("def_structures/def_structures_plastic.jl")
    include("def_structures/def_structures_DFT/def_structures_DFT_SEP.jl")
    include("def_structures/def_structures_DFT/def_structures_DFT_EP.jl")
    include("def_structures/def_structures_DFT/def_structures_DFT_PSE.jl")
    include("def_structures/def_structures_DFT/def_structures_DFT_PP.jl")
    include("def_structures/def_structures_DFT/def_structures_DFT_PE.jl")
    include("def_structures/def_structures.jl")
    include("interactions/interaction_beam_to_rigid/def_interactions.jl")
    include("sparse_matrices.jl")

    include("beams/def_node_beam.jl")
    include("beams/def_beam/def_beam_elastic.jl")
    include("beams/def_beam/def_beam_superelastic.jl")
    include("beams/def_beam/def_beam_plastic.jl")
    include("beams/def_beam/def_beam_DFT/def_beam_DFT_EP.jl")
    include("beams/def_beam/def_beam_DFT/def_beam_DFT_PE.jl") 
    include("beams/def_beam/def_beam_DFT/def_beam_DFT_SEP.jl")
    include("beams/def_beam/def_beam_DFT/def_beam_DFT_PSE.jl")
    include("beams/def_beam/def_beam_DFT/def_beam_DFT_PP.jl")
    include("beams/def_beam/def_beam.jl")

    include("beams/compute_beams.jl")
    include("beams/assemble_beams.jl")
    include("beams/utils/utils_beams.jl")
    include("interactions/interaction_beam_to_rigid/utils_compute_interactions_beams.jl")
    include("interactions/interaction_beam_to_rigid/interactions_beams_regularization.jl")
    include("interactions/interaction_beam_to_rigid/compute&assemble_interactions_beams.jl")

    include("beams/utils/utils_tensors.jl")
    include("beams/materials/linear_elastic_beam/linear_elastic_beam.jl")
    include("beams/materials/superelastic_beam/superelastic_beam.jl")
    include("beams/materials/superelastic_beam/utils_superelastic.jl")
    include("beams/materials/platic_beam/utils_plastic.jl")
    include("beams/materials/platic_beam/plastic_beam.jl")
    include("beams/materials/DFT_beam/DFT_EP_beam.jl")
    include("beams/materials/DFT_beam/DFT_SEP_beam.jl")
    include("beams/materials/DFT_beam/DFT_PSE_beam.jl")
    include("beams/materials/DFT_beam/DFT_PP_beam.jl")
    include("beams/materials/DFT_beam/DFT_PE_beam.jl")

    include("inizialization.jl")
    include("run_simulation.jl")
    include("solve_step_dynamics.jl")
    include("predictor.jl")
    include("corrector.jl")
    include("utils_bcs_and_loads.jl")
    include("IO/write_output_files.jl")
    include("utils_solver.jl")

    include("IO/visualization_beams.jl") 
    include("IO/visualization.jl") 
    include("IO/read_input_files.jl")  

    export NodesBeams, BeamsConfiguration
    export ElasticBeams, SuperElasticBeams, DFT_SEPBeams, DFT_EPBeams, DFT_PSEBeams, DFT_PPBeams, DFT_PEBeams, PlasticBeams
    export ConcentratedForce, Loads
    export Encastre, ImposedDisplacement, BoundaryConditions
    export SimulationParams, run_simulation!
    export read_vtk_tetrahedral_mesh, read_vtk_triangle_mesh
    export RigidInteraction, MultiRigidInteraction, SoftInteraction, PlaneSurface, MovingPlaneSurface, SphereSurface, CylinderSurface, MovingCylinderSurface, MovingSphereSurface, DiscreteSignedDistanceField, TriangulatedSurface, BeamElementSurface, BeamNodeSurface, InteractionProperties
    export Matrices
    
    include("precompile/precompiles.jl")
    _precompile_()

end 