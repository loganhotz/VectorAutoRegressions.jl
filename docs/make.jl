using Documenter
using VectorAutoRegressions

makedocs(
    sitename = "VectorAutoRegressions",
    pages = ["Introduction" => "index.md",
             "API"          => "api.md"],
    # pages = [
    #          "Introduction" => "index.md",
    #          "Models" => "model.md", # covers models, formulas, structural
    #          "Impulse Responses" => "impulse.md",
    #          "Utilities" => "utilities.md"
    # ],
    format = Documenter.HTML(),
    modules = [VectorAutoRegressions],
    checkdocs = :exports
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
