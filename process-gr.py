import sys

if __name__ == "__main__":
    args = sys.argv

    print(args)
    input_graph = args[1]
    output_graph = args[2]

    nyfile = open(input_graph, 'r')

    outfile = open(output_graph, 'w+')

    for line in nyfile:
        linef = line.replace('a ', '')
        outfile.write(linef)

    nyfile.close()
    outfile.close()