
import math


def cal_tm_santalucia(temp_gene):
    """
    计算一小段基因(temp_gene)的tm
    :param temp_gene:
    :return: 这段基因的tm
    """
    AATT = ATTA = TAAT = CAGT = GTCA = CTGA = GACT = CGGC = GCCG = GGCC = 0
    for i in range(len(temp_gene) - 1):
        if (temp_gene[i:i + 2] == 'AA') | (temp_gene[i:i + 2] == 'TT'):
            AATT += 1
        elif temp_gene[i:i + 2] == 'AT':
            ATTA += 1
        elif temp_gene[i:i + 2] == 'TA':
            TAAT += 1
        elif (temp_gene[i:i + 2] == 'CA') | (temp_gene[i:i + 2] == 'TG'):
            CAGT += 1
        elif (temp_gene[i:i + 2] == 'GT') | (temp_gene[i:i + 2] == 'AC'):
            GTCA += 1
        elif (temp_gene[i:i + 2] == 'CT') | (temp_gene[i:i + 2] == 'AG'):
            CTGA += 1
        elif (temp_gene[i:i + 2] == 'GA') | (temp_gene[i:i + 2] == 'TC'):
            GACT += 1
        elif temp_gene[i:i + 2] == 'CG':
            CGGC += 1
        elif temp_gene[i:i + 2] == 'GC':
            GCCG += 1
        elif (temp_gene[i:i + 2] == 'GG') | (temp_gene[i:i + 2] == 'CC'):
            GGCC += 1



    H = AATT * (-7.9-8.4) + ATTA * (-7.2-6.5) + TAAT * (-7.2-6.3) + CAGT * (-8.5-7.4) + GTCA * (-8.4-8.6) + CTGA * (-7.8-6.1) + GACT * (
        -8.2-7.7) + CGGC * (-10.6-10.1) + GCCG * (-9.8-11.1) + GGCC * (-8.0-6.7) + 0.1 + 2.3
    S = AATT * (-22.2-23.6) + ATTA * (-20.4-18.8) + TAAT * (-21.3-18.5) + CAGT * (-22.7-19.3) + GTCA * (-22.4-23.0) + CTGA * (
        -21.0-16.1) + GACT * (
            -22.2-20.3) + CGGC * (-27.2-25.5) + GCCG * (-24.4-28.4) + GGCC * (-19.9-15.6) - 2.8 + 4.1 - 1.4 -5.9-9.0-1.4

    c_Na = 1  # mmol /

    c_K = 50 / 1000  #
    c_Mg = 2.2 / 1000  #
    c_dNTPs = 0.2 / 1000
    c_Tris = 10 / 1000  # mol / L#

    c_oligo = 10 / 1e9  # 寡核苷酸
    c_t = 500 / 1e9  # 引物cd

    # tm = (H * 1000) / (S + 1.987 * math.log((c_t / 1000) / 4)) + 16.6 * math.log(c_Na) - 273.15
    # print(tm)
    c_Mon = c_K + c_Tris  # + c_Na
    c_Mg = c_Mg - c_dNTPs

    kelvins = 273.15

    a = 3.92e-5
    b = -9.11e-6
    c = 6.26e-5
    d = 1.42e-5
    e = -4.82e-4
    f = 5.25e-4
    g = 8.31e-5

    n_bp = len(temp_gene)
    f_GC = (temp_gene.count("C") + temp_gene.count("G")) / n_bp  # 计算一小片段中gc的含量

    tm = (H * 1000) / (S + 1.987 * math.log((c_t / 1000) / 4)) + 16.6 * math.log(c_Na)

    if c_Mon == 0:
        tem = 1 / tm + a + b * math.log(c_Mg) + f_GC * (c + d * math.log(c_Mg)) + (
                e + f * math.log(c_Mg) + g * (math.log(c_Mg) ** 2)) / (2 * (n_bp - 1))
        return 1 / tem - kelvins
    else:
        R = math.sqrt(c_Mg) / c_Mon

        if R < 0.22:
            c_Mon = c_K + c_Tris + c_Na
            tem = 1 / tm + (4.29 * f_GC - 3.95) * 10e-5 * math.log(c_Mon) + 9.4e-6 * (math.log(c_Mon)) ** 2
            return 1 / tem - kelvins

        elif R < 6.0:
            a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(c_Mon) * math.log(c_Mon))
            d = 1.42e-5 * (1.279 - 4.03e-3 * math.log(c_Mon) - 8.03e-3 * (math.log(c_Mon)) ** 2)
            g = 8.31e-5 * (0.486 - 0.258 * math.log(c_Mon) + 5.25e-3 * (math.log(c_Mon)) ** 3)

            tem = 1 / tm + a + b * math.log(c_Mg) + f_GC * (c + d * math.log(c_Mg)) + (
                    e + f * math.log(c_Mg) + g * (math.log(c_Mg) ** 2)) / (2 * (n_bp - 1))
            return 1 / tem - kelvins

        else:
            tem = 1 / tm + a + b * math.log(c_Mg) + f_GC * (c + d * math.log(c_Mg)) + (
                    e + f * math.log(c_Mg) + g * (math.log(c_Mg) ** 2)) / (2 * (n_bp - 1))
            return 1 / tem - kelvins



if __name__ == '__main__':
    str1 = "TAAGCACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTATAGTCGAATAA"[-27:]
    str2 = "GGTAAAATAGTCAACACGCACGGTGTTATTCGACTATAACAAACCATTTTCTTGCGT"[:32]

    f = 'TAAGCACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTATAGTCGAATAA'
    r = 'GGTAAAATAGTCAACACGCACGGTGTTATTCGACTATAACAAACCATTTTCTTGCGT'[::-1]
    f1 = 'CACCGTGCGTGTTGACTATTTTACCTCTGGCGGTGATATACTAGAGAAAGAGG'

    print(f)
    print(' '*26, r)
    print((len(f)-1)*' ', f1)
    print(len('GTGGCACGCACAACTGATAAAATGG'))

    print(str1, cal_tm_santalucia(str1))
    print()
    print(str2, cal_tm_santalucia(str2))