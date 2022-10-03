
# Thomas Freeman 2/11/22
# implementation of the Needleman-Wunsch algorithm
# to count gaps in strings.
import logging
import argparse


def parse(fname):
    with open(fname) as fp:
        lines = fp.readlines()

    return lines[0].rstrip(), lines[1]


def get_args():
    parser = argparse.ArgumentParser(description='Global sequence alignment')

    parser.add_argument('-i', '--input', dest='filename', type=str, help='Input filename', default="./test.txt")
    parser.add_argument('-a', '--algorithm', dest='algorithm', type=str, default="global")
    parser.add_argument('-m', '--match', dest='match_score', type=int, help='Match score', default=20)
    parser.add_argument('-p', '--miss', dest='miss_penalty', type=int, help='Miss match penalty', default=-10)
    parser.add_argument('-o', '--gap_opening', dest='gap_opening', type=int, help='Gap opening penalty', default=-5)
    parser.add_argument('-e', '--gap_extension', dest='gap_extension', type=int, help='Gap extension penalty', default=-5)

    args = parser.parse_args()
    return args

def init_matrix(rows, columns, algorithm='global'):
    mat = [[0] * columns for row in range(rows)]
    gap_mat = [[[0, set()]] * columns for row in range(rows)]
    gap_mat[1][0] = [0, set()]
    gap_mat[0][1] = [0, set()]

    # fills the cells in the first row and first column 
    for i in range(1, rows):
        mat[i][0] = mat[i - 1][0] + gap_penalty(gap_mat, i - 1, 0, 'v')
        gap_mat[i][0] = [1, {'v'}]
    for j in range(1, columns):
        mat[0][j] = mat[0][j - 1] + gap_penalty(gap_mat, 0, j - 1, 'h')
        gap_mat[0][j] = [1, {'h'}]
    
    return mat, gap_mat

def print_matrix(matrix):
    for row in matrix:
        print(''.join(['{0:>{w}}'.format(item, w=5) for item in row]), end='\n\n')

def is_match(char_a, char_b):
    return match_score if char_a == char_b else miss_match_penalty

def gap_penalty(gap_matrix, row_index, col_index, gap_direction):
    
    # If there's a gap detected, we determine if
    # this gap is between both strings or
    # between characters of one string.
    if gap_matrix[row_index][col_index][0] == 1:
        if gap_direction == 'v':
            if 'v' in gap_matrix[row_index][col_index][1]:
                return gap_extend
        elif gap_direction == 'h':
            if 'h' in gap_matrix[row_index][col_index][1]:
                return gap_extend

    return gap_opening

    # Method for the Needleman-Wunsch Algorithm
def global_alignment(M, gap_matrix, a, b, rows, columns):
    
    # fills the score matrix
    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = M[i - 1][j - 1] + is_match(a[j - 1], b[i - 1])
            vgap = M[i - 1][j] + gap_penalty(gap_matrix, i - 1, j, 'v')
            hgap = M[i][j - 1] + gap_penalty(gap_matrix, i, j - 1, 'h')

            options = [diagonal, vgap, hgap]
            index_max = options.index(max(options))

            if options[index_max] == hgap:
                gap_matrix[i][j] = [1, {'h'}]

            if options[index_max] == vgap:
                gap_matrix[i][j] = [1, {'v'}]

            M[i][j] = options[index_max]

    i, j = rows - 1, columns - 1
    aligned_a, aligned_b, mid = (' ') * 3

    # backtrack process from lower right to upper left of the matrix
    while i > 0 and j > 0:
        diagonal = M[i][j] - is_match(a[j - 1], b[i - 1])
        vgap = M[i][j] - gap_penalty(gap_matrix, i - 1, j, 'v')
        hgap = M[i][j] - gap_penalty(gap_matrix, i, j - 1, 'h')

        if M[i - 1][j - 1] == diagonal:
            aligned_a += a[j - 1]
            aligned_b += b[i - 1]
            if is_match(a[j - 1], b[i - 1]) == match_score:
                mid += '|'
            else:
                mid += ' '
            i = i - 1
            j = j - 1
        elif M[i - 1][j] == vgap:
            aligned_a += '-'
            aligned_b += b[i - 1]
            mid += ' '
            i = i - 1
        elif M[i][j - 1] == hgap:
            aligned_a += a[j - 1]
            aligned_b += '-'
            mid += ' '
            j = j - 1

    while j > 0:
        aligned_a += a[j - 1]
        aligned_b += '-'
        mid += ' '
        j = j - 1

    while i > 0:
        aligned_a += '-'
        aligned_b += b[i - 1]
        mid += ' '
        i = i - 1
        
    # The Aligned arrays are printed.
    print(aligned_a[::-1] + '\n' + mid[::-1] + '\n' + aligned_b[::-1], "\n")

    match_hit = mid.count('|')

    return M, M[rows - 1][columns - 1]

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    args = get_args()

    if args.filename:

        filename = args.filename
        sequence_a, sequence_b = parse(filename)

        rows, columns = len(sequence_b) + 1, len(sequence_a) + 1

        match_score = args.match_score
        miss_match_penalty = args.miss_penalty
        gap_opening = args.gap_opening
        gap_extend = args.gap_extension

        print("Filename\t:{}\nAlgorithm\t:{}\nMatch score\t:{}\nMiss match\t:{}\nGap opening\t:{}\nGap Extension\t:{}\n".format
              (args.filename, args.algorithm, match_score, miss_match_penalty, gap_opening, gap_extend))

        if args.algorithm == 'global':
            D, gap_matrix = init_matrix(rows, columns, algorithm=args.algorithm)
            D, score = global_alignment(D, gap_matrix, sequence_a, sequence_b, rows, columns)
            print("Total Alignment Score:", score)
        
