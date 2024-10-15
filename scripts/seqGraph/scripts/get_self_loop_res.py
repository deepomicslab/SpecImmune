import sys

if __name__ == "__main__":
    res = sys.argv[1]
    self_loop_file = open(sys.argv[2], 'w')

    with open(res) as all_f:
        while True:
            line1 = all_f.readline().strip()
            if not line1:
                break
            if line1.startswith("self"):
                line2 = all_f.readline().strip()
                if not line2:
                    break
                self_loop_file.write("self loop: \n")
                self_loop_file.write(line2 + "\n")

    self_loop_file.close()
