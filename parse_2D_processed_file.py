# import pathlib

# step_to_parse = "37000"
# files_to_parse = list(pathlib.Path('.').glob(f'*{step_to_parse}*.txt'))

CONST_X = 0
CONST_Y = 10

def parse_file(filename: str = ""):
    file_to_parse = str(filename)
    f = open(f"{file_to_parse[:file_to_parse.index('.txt')][-1]}_processed.txt", "a")
    with open(filename) as file_in:
        for line in file_in:
            splitted = line.replace("\n", "").split(" ")
            values = list(filter(None, splitted))  # i, j, value
            i, j, value = values

            if int(i) == CONST_X:
                f.write("{} {} \n".format(j, value))

    f.close()


if __name__ == '__main__':
    parse_file(filename="u.txt")
