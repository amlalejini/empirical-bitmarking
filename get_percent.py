import csv

filename = "performance.csv"

and_data = and_avg = {
    "naked": {
        "bitset": [],
        "bitvector": []
    },
    "vector": {
        "bitset": [],
        "bitvector": []
    }
}

with open(filename) as file:
    csv_reader = csv.reader(file)
    container = ""
    treatment = ""

    for row in csv_reader:
        if row[3] == "AND":
            if row[-1] == "T":
                container = "naked"
            else:
                container = "vector"

            if row[1] == "bit_set":
                treatment = "bitset"
            else:
                treatment = "bitvector"
            and_data[container][treatment].append(int(row[4]))

for upper in and_avg.keys():
    for key in and_avg[upper].keys():
        and_l = and_data[upper][key]
        and_avg[upper][key] = sum(and_l) / len (and_l)
print(and_avg)

print("bitset/bitvector relative time")

t_per = and_avg["naked"]["bitvector"] / and_avg["naked"]["bitset"]
vect_per = and_avg["vector"]["bitvector"] / and_avg["vector"]["bitset"]

print("'T'        : bitvector takes " + str(round(t_per, 2)) + "x longer")
print("'vector<t>': bitvector takes " + str(round(vect_per, 2)) + "x longer")

print()
print("bitvector speedup")

speedup = and_avg["vector"]["bitvector"] / and_avg["naked"]["bitvector"]

print("it is " + str(round(speedup, 2)) + "x slower to use a vector of bitvectors")