using CSV, DataFrames, Statistics, Bootstrap

# processing data for feeding comparison
cd(joinpath(ENV["HOME"], "Dropbox", "KillifishFeederPaper_AndrewMcKay","Revision", "Code", "GitHub","Fig2"))
raw = DataFrame(CSV.File("Input/Figure2-SourceData3_feedingdata_auto_manual.csv"))
# rename for easier working
raw = rename(raw, "Scaled_to_one_manual_feeding_at_target_20mg" => "ScaledFeed")

dfperson = raw[raw[:Type] .== "Person",:]
dfpreciseperson = raw[raw[:Accuracy] .== "PrecisePerson",:]
dffeeder = raw[raw[:Type] .== "Feeder",:]
dffeeder1 = raw[raw[:Feeder] .== "Feeder4",:]
# seems like dfmean should be here somewhere?
dfmeans = by(raw, :Accuracy, :ScaledFeed => mean)

# fraction deviation
dfperson[:deviation] = (dfperson[:ScaledFeed] .- dfmeans[dfmeans[:Accuracy] .== "Person", :ScaledFeed_mean]) ./ dfmeans[dfmeans[:Accuracy] .== "Person", :ScaledFeed_mean]

dfpreciseperson[:deviation] = (dfpreciseperson[:ScaledFeed] .- dfmeans[dfmeans[:Accuracy] .== "PrecisePerson", :ScaledFeed_mean]) ./ dfmeans[dfmeans[:Accuracy] .== "PrecisePerson", :ScaledFeed_mean]

dffeeder[:deviation] = (dffeeder[:ScaledFeed] .- dfmeans[dfmeans[:Accuracy] .== "Feeder", :ScaledFeed_mean]) ./ dfmeans[dfmeans[:Accuracy] .== "Feeder", :ScaledFeed_mean]

dffeeder1[:deviation] = (dffeeder1[:ScaledFeed] .- dfmeans[dfmeans[:Accuracy] .== "Feeder", :ScaledFeed_mean]) ./  dfmeans[dfmeans[:Accuracy] .== "Feeder", :ScaledFeed_mean]

dfperson[:category] = "Multiple People"
dfpreciseperson[:category] = "Single Person"
dffeeder[:category] = "Automatic Feeder"
dffeeder1[:category] = "Single Automatic Feeder"

# this has more rows than original (duplicates) but is intentional because will only plot later
raw = vcat(dfperson, dffeeder, dfpreciseperson, dffeeder1)

# plot out deviations
CSV.write("Input/rawdeviations_percent.csv", raw)


# creating dataframe for moving to R
n_boot = 10000
df = by(raw, [:category], p -> bootstrap(std, Array(p.ScaledFeed), BasicSampling(n_boot)))
df[:std_est] = map(p -> p.t0[1], df[:x1])
df[:std_error] = map(p -> stderror(p)[1], df[:x1])
df[:count] = map(p -> size(p.data)[1], df[:x1])
df[:precision] = map(p -> (1/(p^2)), df[:std_est])
df[:conf] = map(p -> confint(p[:x1], BasicConfInt(0.68))[1], eachrow(df))
df[:precision_conf] = map(p -> map(pp -> (1/(pp^2)), p), df[:conf])

# since we can't transport bootstrap sampling type over to R
select!(df, Not(:x1))
dfout = df

CSV.write("Output/FeedingReproducibility2.20.19data5.26.20.csv", dfout)
df2 = copy(dfout)
df2 = df[:, [:category]]
x = [[x[1] x[2] x[3]] for x in dfout[:precision_conf]]
x = vcat(x...)
df2[:precision] = x[:,1]
df2[:upper] = x[:,2]
df2[:lower] = x[:,3]
CSV.write("Output/precision_confint.csv", df2)
