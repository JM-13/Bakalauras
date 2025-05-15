import pandas as pd

# line = ' ! R1    R(1,2)                  1.5525         -DE/DX =    0.0                 !'
#
# Table = []
#
# raw_row = line.split()
# row = raw_row[1:4] + [raw_row[6]]
# print(row)
# for val in range(1, len(row)):
#     row[val] = float(row[val])
#
# Table.append(row)
# Table.append(row)
# Table.append(row)
# Table.append(row)
# Table.append(row)
# Table.append(row)
# Table.append(row)
# Table.append(row)
# columns = ["ID", "Value1", "Value2", "Value3", "Value4", "Value5", "Value6"]
# df = pd.DataFrame(Table, columns=columns)
# print(df)

line = 'D(20,1,11,18)'
print(list(line))
