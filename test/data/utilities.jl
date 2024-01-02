DATA_DIR = @__DIR__


function load_data(dataset::Symbol = :fred)
    if dataset === :fred
        return load_fred_data()
    end
end



"""
    load_fred_data()
"""
load_fred_data() = load_fred_data(DATA_DIR)
function load_fred_data(dir_path)
    # read in the quarterly data from FRED
    file_name = "fredgraph.csv"
    file_path = joinpath(DATA_DIR, file_name)
    df        = CSV.read(file_path, DataFrame, missingstring=["."])

    # rename, and then comput inflation rate
    rename!(df, [:date, :gdp_pi, :urate, :ffr])
    df[!, :infl] = 400 * log.( df.gdp_pi ./ lag(df.gdp_pi) )

    # restrict the dates
    first, final = Date("1960-01-01"), Date("2000-12-31")
    idx          = (df.date .≥ first) .& (df.date .≤ final)
    df           = df[idx, :]

    df = df[:, [:date, :infl, :urate, :ffr]]
    return df
end
