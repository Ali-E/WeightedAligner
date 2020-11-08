import sys
import pandas as pd


def reverse(pattern):
    patt_rev_comp = ""
    for idx in range(len(pattern) - 1, -1, -1):
        patt_rev_comp += pattern[idx]
    return patt_rev_comp


def reverse_complement(pattern):
    comp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    patt_rev_comp = ""
    for idx in range(len(pattern) - 1, -1, -1):
        patt_rev_comp += comp_dict[pattern[idx]]
    return patt_rev_comp


def build_scoring_matrix(alphabet, diag_score=1, off_diag_score=-1, dash_score=-1, varying_scores=True, AG=-0.5, CT=-0.75, score_dict={}, reverse=False):
    """
    The function returns a dictionary of dictionaries as the scoring matrix.
    """
    scoring_matrix = {}
    alph_list = list(alphabet) + ['-']
    for alph in alph_list:
        scoring_matrix[alph] = {}
        for alph_2 in alph_list:
            if alph == '-' or alph_2 == '-':
                scoring_matrix[alph][alph_2] = dash_score
            elif alph == alph_2:
                scoring_matrix[alph][alph_2] = diag_score
            else:
                # scoring_matrix[alph][alph_2] = off_diag_score

                if alph == 'N' or alph_2 == 'N':
                    scoring_matrix[alph][alph_2] = 0
                else:
                    scoring_matrix[alph][alph_2] = off_diag_score

                if varying_scores:
                    if reverse:
                        tup = (alph_2, alph)
                    else:
                        tup = (alph, alph_2)
                    if tup[0] == 'G' and tup[1] == 'A':
                        scoring_matrix[alph][alph_2] = AG
                    if tup[0] == 'T' and tup[1] == 'C':
                        scoring_matrix[alph][alph_2] = CT
                    if (tup[1], tup[0]) in score_dict:
                        scoring_matrix[alph][alph_2] = score_dict[(tup[1], tup[0])]

    return scoring_matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    a_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        temp = a_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]['-']
        if global_flag == 0 and temp < 0:
            temp = 0
        a_matrix[idx_x][0] = temp

    for idx_y in range(1, len(seq_y) + 1):
        temp = a_matrix[0][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
        if global_flag == 0 and temp < 0:
            temp = 0
        a_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(a_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-'],
                       a_matrix[idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]],
                       a_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])
            if global_flag == 0 and temp < 0:
                temp = 0
            a_matrix[idx_1][idx_2] = temp

    return a_matrix


def compute_alignment_matrix_1(seq_x, seq_y, scoring_matrix, global_flag):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    computes the max_amount of the matrix.
    """
    max_amount = -1
    a_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        temp = a_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]['-']
        if global_flag == 0 and temp < 0:
            temp = 0
        if temp > max_amount:
            max_amount = temp
        a_matrix[idx_x][0] = temp

    for idx_y in range(1, len(seq_y) + 1):
        temp = a_matrix[0][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
        if global_flag == 0 and temp < 0:
            temp = 0
        if temp > max_amount:
            max_amount = temp
        a_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(a_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-'],
                       a_matrix[idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]],
                       a_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])
            if global_flag == 0 and temp < 0:
                temp = 0
            if temp > max_amount:
                max_amount = temp
            a_matrix[idx_1][idx_2] = temp

    return max_amount


def compute_affine_alignment_matrix(seq_x, seq_y, scoring_matrix, gop, gep, global_flag):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    m_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]
    l_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]
    u_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        temp = max(l_matrix[idx_x - 1][0] - gep, m_matrix[idx_x - 1][0] - gop)
        l_matrix[idx_x][0] = temp

        temp = 0
        u_matrix[idx_x][0] = temp

        temp = max(u_matrix[idx_x][0], l_matrix[idx_x][0])
        m_matrix[idx_x][0] = temp

    for idx_y in range(1, len(seq_y) + 1):
        temp = 0
        l_matrix[0][idx_y] = temp

        temp = max(u_matrix[0][idx_y - 1] - gep, m_matrix[0][idx_y - 1] - gop)
        u_matrix[0][idx_y] = temp

        temp = max(u_matrix[0][idx_y], l_matrix[0][idx_y])
        m_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(l_matrix[idx_1 - 1][idx_2] - gep, m_matrix[idx_1 - 1][idx_2] - gop)
            l_matrix[idx_1][idx_2] = temp

            temp = max(u_matrix[idx_1][idx_2 - 1] - gep, m_matrix[idx_1][idx_2 - 1] - gop)
            u_matrix[idx_1][idx_2] = temp

            temp = max(l_matrix[idx_1][idx_2], u_matrix[idx_1][idx_2],
                       m_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])
            m_matrix[idx_1][idx_2] = temp

    return (l_matrix, m_matrix, u_matrix)


def compute_kg_alignment_matrix(seq_x, seq_y, scoring_matrix, Klimit):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    u_matrix = [[[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)] for dummy_idx_3 in range(Klimit + 1)]
    m_matrix = [[[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)] for dummy_idx_3 in range(Klimit + 1)]
    l_matrix = [[[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)] for dummy_idx_3 in range(Klimit + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        for idx_y in range(0, len(seq_y) + 1):
            temp = l_matrix[0][idx_x - 1][idx_y] + scoring_matrix[seq_x[idx_x - 1]]['-']
            l_matrix[0][idx_x][idx_y] = temp

    for idx_y in range(1, len(seq_y) + 1):
        for idx_x in range(0, len(seq_x) + 1):
            temp = u_matrix[0][idx_x][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
            # print(temp)
            u_matrix[0][idx_x][idx_y] = temp

    for idx_x in range(0, len(seq_x) + 1):
        for idx_y in range(0, len(seq_y) + 1):
            temp = max(u_matrix[0][idx_x][idx_y], l_matrix[0][idx_x][idx_y])
            m_matrix[0][idx_x][idx_y] = temp

    for k in range(1, Klimit+1):
        for idx_x in range(1, len(seq_x) + 1):
            for idx_y in range(0, len(seq_y) + 1):
                temp = max(l_matrix[k][idx_x - 1][idx_y], m_matrix[k-1][idx_x - 1][idx_y]) + scoring_matrix[seq_x[idx_x - 1]]['-']
                l_matrix[k][idx_x][idx_y] = temp

        for idx_y in range(1, len(seq_y) + 1):
            for idx_x in range(0, len(seq_x) + 1):
                temp = max(u_matrix[k][idx_x][idx_y - 1], m_matrix[k-1][idx_x][idx_y-1]) + scoring_matrix['-'][seq_y[idx_y - 1]]
                u_matrix[k][idx_x][idx_y] = temp

        for idx_x in range(0, len(seq_x) + 1):
            for idx_y in range(0, len(seq_y) + 1):
                temp = max(u_matrix[k][idx_x][idx_y], l_matrix[k][idx_x][idx_y], m_matrix[k][idx_x-1][idx_y-1]+scoring_matrix[seq_x[idx_x-1]][seq_y[idx_y-1]])
                m_matrix[k][idx_x][idx_y] = temp

    return l_matrix, m_matrix, u_matrix


def compute_kg_alignment(seq_x, seq_y, Klimit, scoring_matrix, alignment_matrix):
    """
    This function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix with affine penalty.
    """
    state = 1
    l_matrix = alignment_matrix[0]
    m_matrix = alignment_matrix[1]
    u_matrix = alignment_matrix[2]
    x_prim = ''
    y_prim = ''
    idx_1 = len(seq_x)
    idx_2 = len(seq_y)

    max_k = -1
    max_elem = -10000
    for k in range(Klimit+1):
        if m_matrix[k][idx_1][idx_2] > max_elem:
            max_k = k
            max_elem = m_matrix[k][idx_1][idx_2]

    curr_k = max_k
    while idx_1 != 0 and idx_2 != 0:
        if state == 1:
            if m_matrix[curr_k][idx_1][idx_2] == m_matrix[curr_k][idx_1-1][idx_2-1]+scoring_matrix[seq_x[idx_1-1]][seq_y[idx_2-1]]:
                print("yeah")
                # print m_matrix[idx_1][idx_2]
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_1 -= 1
                idx_2 -= 1
            elif m_matrix[curr_k][idx_1][idx_2] == l_matrix[curr_k][idx_1][idx_2]:
                state = 0
                print(state)
            elif m_matrix[curr_k][idx_1][idx_2] == u_matrix[curr_k][idx_1][idx_2]:
                state = 2
                print(state)

            else:
                print("fault occured!")
                return 0
            continue
        elif state == 0:
            if l_matrix[curr_k][idx_1][idx_2] == m_matrix[curr_k-1][idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
                state = 1
            elif l_matrix[curr_k][idx_1][idx_2] == l_matrix[curr_k][idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1

            else:
                print("fault!")
                return 0
            continue
        elif state == 2:

            if u_matrix[curr_k][idx_1][idx_2] == m_matrix[curr_k-1][idx_1][idx_2-1] + scoring_matrix['-'][seq_y[idx_2 - 1]]:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1
                state = 1
            elif u_matrix[curr_k][idx_1][idx_2] == u_matrix[curr_k][idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]]:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1
            else:
                print("fault")
                return 0
            continue

    print(x_prim)
    print(y_prim)
    while idx_1 != 0:
        x_prim = seq_x[idx_1 - 1] + x_prim
        y_prim = '-' + y_prim
        idx_1 -= 1
    while idx_2 != 0:
        x_prim = '-' + x_prim
        y_prim = seq_y[idx_2 - 1] + y_prim
        idx_2 -= 1

    print(max_k)
    return m_matrix[max_k][len(seq_x)][len(seq_y)], x_prim, y_prim


def compute_fit_alignment_matrix(seq_x, seq_y, scoring_matrix):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    a_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        """
        temp = a_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]['-']
        if global_flag == 0 and temp < 0:
            temp = 0
        """
        a_matrix[idx_x][0] = 0

    for idx_y in range(1, len(seq_y) + 1):
        temp = a_matrix[0][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
        a_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(a_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-'],
                       a_matrix[idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]],
                       a_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])

            a_matrix[idx_1][idx_2] = temp

    return a_matrix


def compute_kgaps_alignment_matrix(seq_x, seq_y, scoring_matrix, Klimit):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    U_matrix = [[[0 for dummy_idx_3 in range(Klimit + 1)] for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]
    M_matrix = [[[0 for dummy_idx_3 in range(Klimit + 1)] for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]
    L_matrix = [[[0 for dummy_idx_3 in range(Klimit + 1)] for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        """
        temp = a_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]['-']
        if global_flag == 0 and temp < 0:
            temp = 0
        """
        a_matrix[idx_x][0] = 0

    for idx_y in range(1, len(seq_y) + 1):
        temp = a_matrix[0][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
        a_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(a_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-'],
                       a_matrix[idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]],
                       a_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])

            a_matrix[idx_1][idx_2] = temp

    return a_matrix


def compute_overlap_alignment_matrix(seq_x, seq_y, scoring_matrix):
    """
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    """
    a_matrix = [[0 for dummy_idx_2 in range(len(seq_y) + 1)] for dummy_idx in range(len(seq_x) + 1)]

    for idx_x in range(1, len(seq_x) + 1):
        """
        temp = a_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]['-']
        if global_flag == 0 and temp < 0:
            temp = 0
        """
        a_matrix[idx_x][0] = 0

    for idx_y in range(1, len(seq_y) + 1):
        temp = a_matrix[0][idx_y - 1] + scoring_matrix['-'][seq_y[idx_y - 1]]
        a_matrix[0][idx_y] = temp

    for idx_1 in range(1, len(seq_x) + 1):
        for idx_2 in range(1, len(seq_y) + 1):
            temp = max(a_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-'],
                       a_matrix[idx_1][idx_2 - 1] + scoring_matrix['-'][seq_y[idx_2 - 1]],
                       a_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]])

            a_matrix[idx_1][idx_2] = temp

    return a_matrix


"""
print(str(compute_alignment_matrix('AGC', 'ACT', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
                                              '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
                                              'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, False)))
"""


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.
    """
    x_prim = ''
    y_prim = ''
    idx_1 = len(seq_x)
    idx_2 = len(seq_y)
    while idx_1 != 0 and idx_2 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][
            seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]][
                '-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    while idx_1 != 0:
        x_prim = seq_x[idx_1 - 1] + x_prim
        y_prim = '-' + y_prim
        idx_1 -= 1
    while idx_2 != 0:
        x_prim = '-' + x_prim
        y_prim = seq_y[idx_2 - 1] + y_prim
        idx_2 -= 1

    return (alignment_matrix[len(seq_x)][len(seq_y)], x_prim, y_prim)


"""
scoring_matrix = {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
                                               '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
                                               'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}
align = compute_alignment_matrix('ACGT', 'AGAGT', scoring_matrix, True)
print(align)

print(str( compute_global_alignment('ACGT', 'AGAGT', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
                                               '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
                                               'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, align) ))
"""


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    """
    max_idx_1 = -1
    max_idx_2 = -1
    max_amount = -1
    for id_1 in range(len(alignment_matrix)):
        for id_2 in range(len(alignment_matrix[0])):
            if alignment_matrix[id_1][id_2] > max_amount:
                max_amount = alignment_matrix[id_1][id_2]
                max_idx_1 = id_1
                max_idx_2 = id_2

    # print(max_idx_1)
    # print(max_idx_2)

    x_prim = ''
    y_prim = ''
    idx_1 = max_idx_1
    idx_2 = max_idx_2
    while alignment_matrix[idx_1][idx_2] != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    """
    while idx_1 != 0:
        x_prim = seq_x[idx_1 -1] + x_prim
        y_prim = '-' + y_prim
        idx_1 -= 1
    while idx_2 != 0:
        x_prim = '-' + x_prim
        y_prim = seq_y[idx_2 -1] + y_prim
        idx_2 -= 1
    """

    return (max_amount, x_prim, y_prim)


def compute_fit_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    """
    max_idx_1 = -1
    max_idx_2 = -1
    max_amount = 0 - len(seq_x) - len(seq_y)
    second_max_idx_1 = -1
    second_max_idx_2 = -1
    second_max_amount = 0 - len(seq_x) - len(seq_y)
    third_max_idx_1 = -1
    third_max_idx_2 = -1
    third_max_amount = 0 - len(seq_x) - len(seq_y)
    fourth_max_idx_1 = -1
    fourth_max_idx_2 = -1
    fourth_max_amount = 0 - len(seq_x) - len(seq_y)
    fifth_max_idx_1 = -1
    fifth_max_idx_2 = -1
    fifth_max_amount = 0 - len(seq_x) - len(seq_y)
    sixth_max_idx_1 = -1
    sixth_max_idx_2 = -1
    sixth_max_amount = 0 - len(seq_x) - len(seq_y)
    seventh_max_idx_1 = -1
    seventh_max_idx_2 = -1
    seventh_max_amount = 0 - len(seq_x) - len(seq_y)
    eighth_max_idx_1 = -1
    eighth_max_idx_2 = -1
    eighth_max_amount = 0 - len(seq_x) - len(seq_y)
    nineth_max_idx_1 = -1
    nineth_max_idx_2 = -1
    nineth_max_amount = 0 - len(seq_x) - len(seq_y)
    tenth_max_idx_1 = -1
    tenth_max_idx_2 = -1
    tenth_max_amount = 0 - len(seq_x) - len(seq_y)
    for id_1 in range(len(alignment_matrix)):
        # for id_2 in range(len(alignment_matrix[0])):
        id_2 = len(alignment_matrix[0]) - 1
        if alignment_matrix[id_1][id_2] >= max_amount:
            tenth_max_amount = nineth_max_amount
            tenth_max_idx_1 = nineth_max_idx_1
            tenth_max_idx_2 = nineth_max_idx_2
            nineth_max_amount = eighth_max_amount
            nineth_max_idx_1 = eighth_max_idx_1
            nineth_max_idx_2 = eighth_max_idx_2
            eighth_max_amount = seventh_max_amount
            eighth_max_idx_1 = seventh_max_idx_1
            eighth_max_idx_2 = seventh_max_idx_2
            seventh_max_amount = sixth_max_amount
            seventh_max_idx_1 = sixth_max_idx_1
            seventh_max_idx_2 = sixth_max_idx_2
            sixth_max_amount = fifth_max_amount
            sixth_max_idx_1 = fifth_max_idx_1
            sixth_max_idx_2 = fifth_max_idx_2
            fifth_max_amount = fourth_max_amount
            fifth_max_idx_1 = fourth_max_idx_1
            fifth_max_idx_2 = fourth_max_idx_2
            fourth_max_amount = third_max_amount
            fourth_max_idx_1 = third_max_idx_1
            fourth_max_idx_2 = third_max_idx_2
            third_max_amount = second_max_amount
            third_max_idx_1 = second_max_idx_1
            third_max_idx_2 = second_max_idx_2
            second_max_amount = max_amount
            second_max_idx_1 = max_idx_1
            second_max_idx_2 = max_idx_2
            max_amount = alignment_matrix[id_1][id_2]
            max_idx_1 = id_1
            max_idx_2 = id_2


    all_tups = []
    
    x_prim = ''
    y_prim = ''
    idx_1 = max_idx_1
    idx_2 = max_idx_2
    # idx_1 = second_max_idx_1
    # idx_2 = second_max_idx_2
    # idx_1 = third_max_idx_1
    # idx_2 = third_max_idx_2
    # idx_1 = fourth_max_idx_1
    # idx_2 = fourth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1
    
    tmp_tup = (max_amount, x_prim, y_prim, max_idx_1, max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = second_max_idx_1
    idx_2 = second_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (second_max_amount, x_prim, y_prim, second_max_idx_1, second_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = third_max_idx_1
    idx_2 = third_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (third_max_amount, x_prim, y_prim, third_max_idx_1, third_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = fourth_max_idx_1
    idx_2 = fourth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (fourth_max_amount, x_prim, y_prim, fourth_max_idx_1, fourth_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = fifth_max_idx_1
    idx_2 = fifth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (fifth_max_amount, x_prim, y_prim, fifth_max_idx_1, fifth_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = sixth_max_idx_1
    idx_2 = sixth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (sixth_max_amount, x_prim, y_prim, sixth_max_idx_1, sixth_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = seventh_max_idx_1
    idx_2 = seventh_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (seventh_max_amount, x_prim, y_prim, seventh_max_idx_1, seventh_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = eighth_max_idx_1
    idx_2 = eighth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (eighth_max_amount, x_prim, y_prim, eighth_max_idx_1, eighth_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = nineth_max_idx_1
    idx_2 = nineth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (nineth_max_amount, x_prim, y_prim, nineth_max_idx_1, nineth_max_idx_2)
    all_tups.append(tmp_tup)


    x_prim = ''
    y_prim = ''
    idx_1 = tenth_max_idx_1
    idx_2 = tenth_max_idx_2
    while idx_2 != 0 and idx_1 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]]['-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    tmp_tup = (tenth_max_amount, x_prim, y_prim, tenth_max_idx_1, tenth_max_idx_2)
    all_tups.append(tmp_tup)


    # print(max_amount, second_max_amount, third_max_amount, fourth_max_amount, fifth_max_amount, sixth_max_amount)
    # return (max_amount, x_prim, y_prim, max_idx_1, max_idx_2)
    return all_tups


def compute_overlap_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    """
    max_idx_1 = -1
    max_idx_2 = -1
    max_amount = -1
    id_1 = len(alignment_matrix) - 1

    for id_2 in range(len(alignment_matrix[0])):
        if alignment_matrix[id_1][id_2] > max_amount:
            max_amount = alignment_matrix[id_1][id_2]
            max_idx_1 = id_1
            max_idx_2 = id_2

    x_prim = ''
    y_prim = ''
    idx_1 = max_idx_1
    idx_2 = max_idx_2
    while idx_2 != 0:
        if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2 - 1] + scoring_matrix[seq_x[idx_1 - 1]][
            seq_y[idx_2 - 1]]:
            x_prim = seq_x[idx_1 - 1] + x_prim
            y_prim = seq_y[idx_2 - 1] + y_prim
            idx_1 -= 1
            idx_2 -= 1
        else:
            if alignment_matrix[idx_1][idx_2] == alignment_matrix[idx_1 - 1][idx_2] + scoring_matrix[seq_x[idx_1 - 1]][
                '-']:
                x_prim = seq_x[idx_1 - 1] + x_prim
                y_prim = '-' + y_prim
                idx_1 -= 1
            else:
                x_prim = '-' + x_prim
                y_prim = seq_y[idx_2 - 1] + y_prim
                idx_2 -= 1

    """
    while idx_1 != 0:
        x_prim = seq_x[idx_1 -1] + x_prim
        y_prim = '-' + y_prim
        idx_1 -= 1
    while idx_2 != 0:
        x_prim = '-' + x_prim
        y_prim = seq_y[idx_2 -1] + y_prim
        idx_2 -= 1
    """

    return (max_amount, x_prim, y_prim)


def parse_fasta(ref_file):
    all_ref_seqs = {}
    ref_seq = ""
    prime_key = ""
    flag_first = True
    with open(ref_file) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] != '>':
                if line[-1] == '\n':
                    ref_seq += line[:-1]
                else:
                    ref_seq += line
            else:
                new_key = line[1:-1]
                if flag_first:
                    flag_first = False
                    prime_key = new_key
                else:
                    trimmed_ref_seq = ref_seq.upper().strip('N')
                    if len(trimmed_ref_seq):
                        all_ref_seqs[prime_key] = trimmed_ref_seq 
                    prime_key = new_key
                    ref_seq = ""

    trimmed_ref_seq = ref_seq.upper().strip('N')
    if len(trimmed_ref_seq):
        all_ref_seqs[prime_key] = trimmed_ref_seq
    return all_ref_seqs


def parse_query_file(seq_file):
    all_seqs = []
    with open(seq_file) as f:
        lines = f.readlines()
        for idx in range(0,len(lines),2):
            if lines[idx+1][-1] == '\n':
                all_seqs.append((lines[idx][:-1], lines[idx+1][:-1]))
            else:
                all_seqs.append((lines[idx][:-1], lines[idx+1]))
    return all_seqs


def align(all_ref_seqs, all_seqs, AG=-0.5, CT=-0.75):
    alphabet = set(['A', 'C', 'G', 'T', 'N', 'R', 'D', 'V', 'W', 'Y', 'S', 'H', 'X', 'F', 'Z', 'U', 'X', 'B', 'E', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q'])
    scoring_matrix = build_scoring_matrix(alphabet, 1, -1, -1, AG=AG, CT=CT)
    scoring_matrix_rev = build_scoring_matrix(alphabet, 1, -1, -1, reverse=True, AG=AG, CT=CT)

    all_results = []
    for idx, seq in enumerate(all_seqs):
        # seen_set = set()
        for key, ref_seq in all_ref_seqs.items():

            # print('ref: ', ref_seq)

            # if idx % 1000 == 0:
            #     print(idx)

            # if seq[0] in seen_set:
            #     continue
            # seen_set.add(seq[0])

            if len(seq[1]) > len(ref_seq):
                rev_flag = 'T'
            else:
                rev_flag = 'F'

            try:
                # scoring_matrix = build_scoring_matrix(alphabet, 1, -1, -1)
                if rev_flag == 'T':
                    align_mat = compute_fit_alignment_matrix(seq[1], ref_seq, scoring_matrix_rev)
                    res_all = compute_fit_alignment(seq[1], ref_seq, scoring_matrix, align_mat)
                    for res in res_all:
                        if float(res[0])/len(ref_seq) > 0.9:
                            print(res)
                        all_results.append((seq[0][1:], key, float(res[0])/len(ref_seq), res[0], len(seq[1]), res[2], res[1], res[4], res[3]))
                else:
                    align_mat = compute_fit_alignment_matrix(ref_seq, seq[1], scoring_matrix)
                    res_all = compute_fit_alignment(ref_seq, seq[1], scoring_matrix, align_mat)
                    for res in res_all:
                        if float(res[0])/len(seq[1]) > 0.9:
                            print(res)
                        all_results.append((seq[0][1:], key, float(res[0])/len(seq[1]), res[0], len(seq[1]), res[1], res[2], res[3], res[4]))
            except Exception as e:
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                print(e, seq[0])
                print(seq[1])


    res_df = pd.DataFrame(all_results, columns=['id', 'ref_id', 'norm_score', 'raw_score', 'length', 'align_ref', 'align_seq', 'ref_end', 'seq_end'])
    return res_df



if __name__ == "__main__":
    
    # filename = sys.argv[1]
    # all_ref_seqs = parse_fasta(filename)
    # print(all_ref_seqs)
    # exit(0)




    # seq_x = 'TGGAGGAGAGCCGCGACCGGCTG'
    # seq_y = 'GCGGGGTCGGCGGCGACGT'

    # print(seq_y)
    # seq_y = reverse(seq_y)
    # seq_y = reverse_complement(seq_y)
    # print(seq_y)

    # # scoring_matrix = build_scoring_matrix(alphabet, 1, -1, -1)
    # print(scoring_matrix)
    # align_mat = compute_alignment_matrix(seq_y, seq_x, scoring_matrix, True)
    # for row in align_mat:
    #     print(row)

    # # result = compute_global_alignment(seq_x, seq_y, scoring_matrix, align_mat)
    # # print(result)
    # res = compute_global_alignment(seq_y, seq_x, scoring_matrix, align_mat)
    # print(res)
   
    # exit(0)
   

    ref_file = sys.argv[1]
    seq_file = sys.argv[2]
    # rev_flag = sys.argv[3]

    all_ref_seqs = parse_fasta(ref_file)

    # print(ref_seq)
    all_seqs = parse_query_file(seq_file)

    # print('last seq: ', all_seqs[-1])
    # print(all_seqs)
    # print(all_ref_seqs)
    

    res_df = align(all_ref_seqs, all_seqs)
    
    res_df.to_csv(seq_file.split(".")[0] + "_" + ref_file.split("/")[1] + "_results.csv", sep='\t', index=False)
    # res_df.to_csv(ref_file.split(".")[0] + "_" + seq_file.split('.')[0] + "_results.csv", sep='\t', index=False)
















