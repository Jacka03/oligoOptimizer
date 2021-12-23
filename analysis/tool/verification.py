import time

from nupack import SetSpec, RawStrand, RawComplex, Strand, Complex, Tube, tube_analysis, \
    Model, complex_analysis, complex_concentrations, Domain, TargetStrand, analysis


class Verification:
    def __init__(self, list_g1, list_g2, len1, temp, conc):
        self.list_g1 = list_g1  # 上面的基因片段
        self.list_g2 = list_g2  # 下面的基因片段
        self.gene_list = list_g1 + list_g2

        self.len_g1 = len(list_g1)  # 上面的基因片段数量
        self.len_g2 = len(list_g2)  # 下面的基因片段数量
        self.len_gene = self.len_g1 + self.len_g2

        self.len1 = len1  # 拼接前的基因片段数

        self.c_gene = conc  # 检验的时候反应的基因序列的浓度

        self.first_check = conc / 10  # 第一次验证时的浓度
        self.second_check = 1e-14  # 第二次验证时的浓度
        self.temp = temp  # 验证 时的温度

    def get_strands_tube_tow(self):
        # 获取试管中只有两条基因片段的所有情况
        tubes = []
        count = 1
        for i in range(self.len_gene):
            for j in range(i, self.len_gene):
                strands = {Strand(self.gene_list[i], name="t{0}:L{1}".format(count, i)): self.c_gene,
                           Strand(self.gene_list[j], name="t{0}:R{1}".format(count, j)): self.c_gene}
                tubes.append(Tube(strands=strands, complexes=SetSpec(max_size=2), name='t{0}'.format(count)))
                count += 1
        return tubes

    def get_tube(self, index1):
        # 找到能与之对应的 index
        if index1 < self.len_g1:
            if self.len1 % 2 == 0:  # 边界
                tem_err_list = [index1, index1 + self.len_g1, index1 + self.len_g1 + 1]
            else:
                tem_err_list = [index1, index1 + self.len_g1]
        else:
            if index1 == self.len_g1 + 1:
                tem_err_list = [index1, index1 - self.len_g1]
            elif index1 == self.len_g1 + self.len_g2 and self.len1 % 2 == 0:
                tem_err_list = [index1, index1 - self.len_g1 - 1]
            else:
                tem_err_list = [index1, index1 - self.len_g1 - 1, index1 - self.len_g1]
        return tem_err_list

    def analysis_two(self):
        my_model = Model(material='dna', celsius=self.temp)  # 温度
        tubes = self.get_strands_tube_tow()  # 得到每个试管中都有两条DNA单链
        tube_results = tube_analysis(tubes=tubes, model=my_model)
        all_conc = {}
        for t in tubes:
            for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
                all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度

        all_conc = sorted(all_conc.items(), key=lambda d: d[1], reverse=True)  # 排序

        error = {}  # 怀疑是错配的
        # 验证
        for k, v in all_conc:
            if k.count("+") == 1 and v > self.first_check:  # 将浓度换成输入的浓度
                # 根据：分割k， 然后根据名字具有顺序关系，然后确定是不是正确的配对
                tem_split = k.split('+')
                t1 = int(tem_split[0].split(':')[1][1:])  # 单链序号
                t2 = int(tem_split[1].split(':')[1][1:-1])
                if t1 > t2:
                    t1, t2 = t2, t1

                if t2 - t1 - self.len_g1 not in [-1, 0]:  # 错配
                    ste = '{0},{1}'.format(t1, t2)
                    if ste in error and v < error[ste]:
                        continue
                    else:
                        error[ste] = v
            elif v < self.first_check:  #
                break

        # 找出那两个的相邻的放到一起，然后反应，看下最后的结果
        # error_end = []  # 经过校验后还是错配的
        # for enu in error:
        #     str = enu.split(',')
        #     tem_err_list = [self.get_tube(int(str[0])), self.get_tube(int(str[1]))]
        #     if self.verification_two(tem_err_list):
        #         error_end.append(enu)

        tem_info = {}
        for key, val in error.items():
            arr = key.split(',')
            arr = [int(i) for i in arr]

            if arr[0] >= self.len_g1:
                str1 = 'R{0}'.format(arr[0] - self.len_g1)
            else:
                str1 = 'F{0}'.format(arr[0])

            if arr[1] >= self.len_g1:
                str2 = 'R{0}'.format(arr[1] - self.len_g1)
            else:
                str2 = 'F{0}'.format(arr[1])
            tem_info[str1 + ',' + str2] = val
        # print(len(tem_info), tem_info)

        # print("目标{0},list1:{1},list2:{2}".format(self.len1, self.len_g1, self.len_g2))
        # print("出错的{0},{1}".format(len(error), error))
        # print("最终检测还是出错的{0}".format(error_end))
        # print('总数{0}'.format(len(all_conc)))
        return tem_info


if __name__ == '__main__':

    data1 = ['AGCACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTTATAGTCGA', 'ATAACACCGTGCGTGTTGACTATTTTACCTCTGGCGGTGATATACTAGAGA', 'AAGAGGAGAAATACTAGATGACCATGATTACGCCAAGCGCGCAATTAAC', 'CCTCACTAAAGGGAACAAAAGCTGGAGCTCCACCGCGGTGGC', 'GGCAGCACTAGAGCTAGTGGATCCCCCGGGCTGTAGAAATTCGATAT', 'CAAGCTTATCGATACCGTCGACCTCGAGGGGGGGCCCGG', 'ACCCAATTCGCCCTATAGTGAGTCGTATTACGCGCGCTCACTGGC', 'CGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCA', 'CTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCG', 'GCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGC', 'GCAGCCTGAATAATAACGCTGATAGTGCTAGTGTAGATCGCTACTAGAGCCAG', 'GCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTG', 'TTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTC', 'GGGTGGGCCTTTCTGCGTT']
    data2 = ['AAAATAGTCAACACGCACGGTGTTATTCGACTATAACAAACCATTTTCTTGCGTAA', 'AATCATGGTCATCTAGTATTTCTCCTCTTTCTCTAGTATATCACCGCCAGAGGT', 'CCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTGCGCGCTTGGCGT', 'GATCCACTAGCTCTAGTGCTGCCGCCACCGCGGTGGAGC', 'AGGTCGACGGTATCGATAAGCTTGATATCGAATTTCTACAGCCCGGGG', 'GACTCACTATAGGGCGAATTGGGTACCGGGCCCCCCCTCG', 'CAGTCACGACGTTGTAAAACGACGGCCAGTGAGCGCGCGTAA', 'GGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCC', 'GGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGG', 'CACTATCAGCGTTATTATTCAGGCTGCGCAACTGTTGGGAAGGGCGA', 'GACTGAGCCTTTCGTTTTATTTGATGCCTGGCTCTAGTAGCGATCTACACTAG', 'AGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTT', 'AACGCAGAAAGGCCCACCCGAAGGTGAGCCAGTGTGACTCTAG']
    data3 = 27
    data4 = 37
    data5 = 1e-8

    veri = Verification(data1, data2, data3, data4, data5)
    ver_info = veri.analysis_two()
    print(len(ver_info), ver_info)


