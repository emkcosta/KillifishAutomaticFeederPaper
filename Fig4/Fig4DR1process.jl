using DataFrames, CSV, Statistics, Dates

function convertdays2weeks(x)
    return Int.(round.(x/7, RoundUp))
end

function convertdays23weeks(x)
    return Int.(round.(x/21, RoundUp))
end

cd(joinpath(ENV["HOME"], "Dropbox", "McKay_CodeChecking", "Code", "Fig4"))

df = CSV.read("../../Data/Table S11.csv")

# make sure all missing observed are unobserved
df[ismissing.(df[:Observed]), :Observed] = 0

# for all without death dates, make today and unobserved
x = (Dates.today())
x = map(p -> p.value, [Month(x) Day(x) Year(x)-Year(2000)])
y = join(x, "/")
df[ismissing.(df[:DeathDate]), :DeathDate] .= y
df[:DeathDate] = map(p -> Dates.Date(p, "m/dd/yy") + Year(2000), df[:DeathDate])

# and calculate lifespan
# put in temp hatch date for some:
df[:HatchDate] = map(p -> Dates.Date(p, "m/dd/yy") + Year(2000), df[:HatchDate])

#df[:DeathDate] = map(p -> Dates.Date(p, "mm/dd/yy") + Year(2000), df[df[:DeathDate], :DeathDate])
df[:Lifespan_days] = df[:DeathDate] .- df[:HatchDate]
df[:Lifespan_days] = map(p -> p.value, df[:Lifespan_days])

#remove Array{Missing,1} from df
df = df[:,map(p -> typeof(df[p]) != Array{Missing,1}, names(df))]

df[:Lifespan_weeks] = convertdays2weeks(df[:Lifespan_days])
df[:Lifespan_3weeks] = convertdays23weeks(df[:Lifespan_days])

# export for lifespan alone
CSV.write("firstlifespandata.csv", df)
