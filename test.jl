using CSV, DataFrames
df = DataFrame(CSV.File("earth_data.csv", delim=','))
