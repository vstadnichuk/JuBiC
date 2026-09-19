include("mibs_aux_reader.jl")
include("mibs_mps_reader.jl")

struct UnsupportedInstanceFormatError <: Exception
    message::String
end

Base.showerror(io::IO, err::UnsupportedInstanceFormatError) = print(io, err.message)

function _validate_supported_instance_format(mps_data::MPSData, aux_data::AUXData)
    leader_rows = vcat(mps_data.rows_less_than, mps_data.rows_greater_than, mps_data.rows_equal)
    for row in leader_rows
        row in aux_data.constraints && continue
        follower_vars = sort!(unique([
            var_name for (var_name, _) in get(mps_data.columns, row, Tuple{String,Number}[])
            if var_name in aux_data.variables
        ]))
        isempty(follower_vars) && continue
        throw(UnsupportedInstanceFormatError(
            "Unsupported BOBILib format: leader constraint '$row' contains follower variable(s) " *
            "$(join(follower_vars, ", ")). JuBiC currently does not support second-level variables in first-level constraints."
        ))
    end
    return nothing
end

include("mibs_instance_builder.jl")
include("mibs_instance_gbc_preprocess.jl")
include("mibs_instance_gbc_builder.jl")
include("mibs_transformation.jl")
include("mibs_instance_to_file_input.jl")

export get_MibS_instance, get_GBC_instance, transform_GBC_to_MibS, output_MibS_instance, output_GBC_instance
