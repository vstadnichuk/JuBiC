import Dates, CSV, DataFrames

function create_folder_if_not_exists(path_to_folder)
    if !isdir(path_to_folder)
        mkpath(path_to_folder)  # Creates the directory and any necessary parent directories.
    end
end


function init_logging_folder()
    # Get the current date and time and format it for the folder name
    current_time = Dates.now()
    formatted_time = Dates.format(current_time, "yyyy-mm-dd_HH-MM-SS")

    # Create the folder name
    folder_name = "logs/$(formatted_time)"

    # Ensure the logs directory does not exists
    if isdir(folder_name)
        error("Folder already exists: $folder_name")
    end

    # Create the directory
    mkpath(folder_name)

    return folder_name
end

"""
    print_stats_to_csv(stats_list, file)

Prints a colllection of RunStats objects into a csv file.
_stats_list_: Collection of RunStats.
_file_: Full path to .csv file. Please do not forget to give the file the correct ending. 
"""
function print_stats_to_csv(stats_list, file)
    df = stats_to_dataframe(stats_list)
    CSV.write(file, df; delim=';', decimal=',')
end


"""
    stats_to_dataframe(stats_list)

Auxilliary function for nice output into csv files of multiple statistics.
"""
function stats_to_dataframe(stats_list)
    all_keys = collect_all_keys(stats_list)
    df = DataFrames.DataFrame()
    for key in all_keys
        df[!, key] = get.([s.data for s in stats_list], key, "")
    end
    return df
end


"""
    collect_all_keys(stats_list)

Auxilliary function for nice output into csv files of multiple statistics.
"""
function collect_all_keys(stats_list)
    all_keys = Set{String}()
    for stats in stats_list
        union!(all_keys, keys(stats.data))
    end
    return sort(collect(all_keys))  # Sort for consistent ordering
end
