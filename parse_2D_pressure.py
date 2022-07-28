import pathlib
from collections import defaultdict

step_to_parse = "30"
files_to_parse = list(pathlib.Path('.').glob(f'*{step_to_parse}*.txt'))

p1 = defaultdict(list)

p2 = defaultdict(list)

def parse_file(filename1, filename2):
    i = [0, 0]
    # j = [_, 0]
    value = 0.0
    k=0
    with open(filename1) as file_in:
        for line in file_in:
            try:
                value = float(line)
                # print(i[1], j[1], k, value)
                # f.write("{} {} {} \n".format(i, k, value))
                k += 1
                p1[i].append(value)

            except Exception as e:
                k = 0
                if "i=" in line:
                    splitted = line.replace("\n", "").split("=")
                    i = splitted[1]
                    # j = splitted[1].replace("\n","").split("=")

    with open(filename2) as file_in:
        for line in file_in:
            try:
                value = float(line)
                # print(i[1], j[1], k, value)
                # f.write("{} {} {} \n".format(i, k, value))
                k += 1
                p2[i].append(value)

            except Exception as e:
                k = 0
                if "i=" in line:
                    splitted = line.replace("\n", "").split("=")
                    i = splitted[1]
                    # j = splitted[1].replace("\n","").split("=")

    # print("")
    # f.close()
    f = open("delta_p.txt", "a")
    for i, v in p1.items():
        f.write("i={}\n".format(i))
        res = [x1 - x2 for (x1, x2) in zip(v, p2[i])]
        for elem in res:
            f.write("{} \n".format(round(elem, 2)))
        f.write("\n\n")




if __name__ == '__main__':
    parse_file(filename1="30-05 13.37 step=100  whole_p.txt", filename2="30-05 14.24 step=200  whole_p.txt")
