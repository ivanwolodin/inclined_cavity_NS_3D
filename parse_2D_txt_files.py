import pathlib

step_to_parse = "37000"
files_to_parse = list(pathlib.Path('.').glob(f'*{step_to_parse}*.txt'))


def parse_file(filename: str = ""):
    file_to_parse = str(filename)
    f = open(f"{file_to_parse[:file_to_parse.index('.txt')][-1]}.txt", "a")
    i = [0, 0]
    # j = [_, 0]
    value = 0.0
    k=0
    with open(filename) as file_in:
        for line in file_in:
            try:
                value = float(line)
                # print(i[1], j[1], k, value)
                f.write("{} {} {} \n".format(i, k, value))
                #                print("{} {} {} \n".format(i[1], k, value))
                k += 1

            except Exception as e:
                k = 0
                if "i=" in line:
                    splitted = line.replace("\n", "").split("=")
                    i = splitted[1]
                    # j = splitted[1].replace("\n","").split("=")

    f.close()


if __name__ == '__main__':
    for item in files_to_parse:
        print (item)
        parse_file(filename=item)
