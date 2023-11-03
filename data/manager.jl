using CSV
using DataFrames
using Dates


DATA_DIR = @__DIR__



function load_data(data::Symbol = :sw2001)
    if data === :sw2001
        return CSV.read(joinpath([DATA_DIR, "stock_watson_2001.csv"]), DataFrame)
    elseif data === :lutkepohl
        return CSV.read(joinpath([DATA_DIR, "lutkepohl.csv"]), DataFrame)
    else
        error("unrecognized dataset code: $data")
    end
end
