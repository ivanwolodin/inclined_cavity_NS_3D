
FILE_TO_PARSE = ""
FILE_TO_WRITE = ""

def parse_file(filename: str = ""):
    f = open(FILE_TO_WRITE, "a")
    i = [_, 0]
    # j = [_, 0]
    value = 0.0
    with open(filename) as file_in:
        for line in file_in:
            try:
                value = float(line)
                # print(i[1], j[1], k, value)
                k += 1
                f.write("{} {} {} \n".format(i[1], k, value))

            except Exception as e:
                k = 0
                if "i=" in line:
                    splitted = line.split("=")
                    i = splitted[0]
                    # j = splitted[1].replace("\n","").split("=")
            finally:
                pass

    f.close()


if __name__ == '__main__':
    parse_file(filename=FILE_TO_PARSE)
