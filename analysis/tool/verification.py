import time

from nupack import SetSpec, RawStrand, RawComplex, Strand, Complex, Tube, tube_analysis, \
    Model, complex_analysis, complex_concentrations, Domain, TargetStrand, analysis


class Verification:

    def __init__(self, oligo, temp, oligo_conc, primer_conc):
        self.oligo = oligo
        # 温度
        self.temp = temp
        # 浓度
        self.oligo_conc = oligo_conc
        self.primer_conc = primer_conc

        self.oligo_check_conc = oligo_conc / 1000  # 第一次验证时的浓度
        self.primer_check_conc = primer_conc / 1000  # 第一次验证时的浓度

    def get_strand_tube_all(self):
        # 获取试管中只有两条基因片段的所有情况

        strands = {}
        for i in range(len(self.oligo)-2):
            strands[Strand(self.oligo[i][1], name=self.oligo[i][0])] = self.oligo_conc

        strands[Strand(self.oligo[-1][1], name=self.oligo[-1][0])] = self.primer_conc
        strands[Strand(self.oligo[-2][1], name=self.oligo[-2][0])] = self.primer_conc

        my_model = Model(material='dna', celsius=self.temp)
        t = Tube(strands=strands, complexes=SetSpec(max_size=2), name='t')  # complexes defaults to [A, B]
        tube_results = tube_analysis(tubes=[t], model=my_model)

        all_conc = {}
        for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
            all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度
        all_conc = sorted(all_conc.items(), key=lambda d: d[1], reverse=True)  # 排序

        # print(all_conc)
        # F last one
        above_oligo_count = int(len(self.oligo)/2) - 2
        tem_char = 'F{0}'.format(above_oligo_count)

        error = {}
        for key, val in all_conc:
            if key.count('+') == 1 and val > self.oligo_check_conc:
                tem_split = key[1:-1].split('+')  # # 先去除括号，在根据+分割字符串

                t1 = tem_split[0][1:]
                t2 = tem_split[1][1:]

                if not t1.isdigit() and not t2.isdigit():
                    if val > self.primer_check_conc:
                        # primer + pirmer
                        error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val

                elif not t1.isdigit():
                    # primer + x, x + primer
                    if tem_split[0][0] == 'R' and tem_split[1] != tem_char:
                        error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                    elif tem_split[0][0] == 'F':
                        error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val

                elif not t2.isdigit():
                    if tem_split[1][0] == 'R' and tem_split[0] != tem_char:
                        error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                    elif tem_split[1][0] == 'F':
                        error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val

                elif tem_split[0][0] == tem_split[1][0]:  # 错配
                    # Rx+Rx ; Fx+Fx
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                elif tem_split[0][0] == 'F' and int(t1) - int(t2) not in [0, 1]:
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                elif tem_split[0][0] == 'R' and int(t1) - int(t2) not in [0, -1]:
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
            # 判断条件有待改进
            elif val < self.oligo_check_conc:  #
                break
        # print(error)
        return error


class Verification1:

    def __init__(self, list_g1, list_g2, len1, temp, conc, primer):
        self.list_g1 = list_g1  # 上面的基因片段
        self.list_g2 = list_g2  # 下面的基因片段
        self.gene_list = list_g1 + list_g2

        self.len_g1 = len(list_g1)  # 上面的基因片段数量
        self.len_g2 = len(list_g2)  # 下面的基因片段数量
        self.len_gene = self.len_g1 + self.len_g2

        self.len1 = len1  # 拼接前的基因片段数

        self.c_gene = conc  # 检验的时候反应的基因序列的浓度

        self.first_check = conc / 1000 # 第一次验证时的浓度
        self.second_check = 1e-14  # 第二次验证时的浓度
        self.temp = temp  # 验证 时的温度

        self.n = 8000  # 每次进行运算的tube个数
        self.F_Primer = primer['F_Primer']
        self.R_Primer = primer['R_Primer']

    def get_strands_tube_tow(self):
        # 获取试管中只有两条基因片段的所有情况
        tubes = []
        count = 1  # 记录tube个数
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
        # start = time.time()
        # tube_results = tube_analysis(tubes=tubes, model=my_model)
        # print("analysis time:{0}".format(time.time() - start))
        all_conc = {}

        tubes_list = [tubes[i:i + self.n] for i in range(0, len(tubes), self.n)]

        for i in tubes_list:
            tube_results = tube_analysis(tubes=i, model=my_model)
            for t in i:
                for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
                    all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度

        # all_conc = {}
        # for t in tubes:
        #     for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
        #         all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度
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

    def get_strands_tube_three(self):
        tubes = []
        count = 1  # 记录试管数目
        for i in range(self.len_gene):
            for j in range(i, self.len_gene):
                for k in range(j, self.len_gene):
                    strands = {Strand(self.gene_list[i], name="t{0}:L{1}".format(count, i)): self.c_gene,
                               Strand(self.gene_list[j], name="t{0}:M{1}".format(count, j)): self.c_gene,
                               Strand(self.gene_list[k], name="t{0}:R{1}".format(count, k)): self.c_gene}
                    tubes.append(Tube(strands=strands, complexes=SetSpec(max_size=3),  # max_size,使用变量代替？
                                      name='t{0}'.format(count)))  # complexes defaults to [A, B]
                    count += 1
        return tubes

    def analysis_three(self):
        my_model = Model(material='dna', celsius=self.temp)
        tubes = self.get_strands_tube_three()  # 得到每个试管中都有两条DNA单链
        start = time.time()
        # 放在一个池时间复杂度太高了
        # tube_results = tube_analysis(tubes=tubes, model=my_model)

        all_conc = {}
        # 将包含上万个tube的list分割成多个包含self.n个tube的二维list
        tubes_list = [tubes[i:i + self.n] for i in range(0, len(tubes), self.n)]
        for i in tubes_list:
            tube_results = tube_analysis(tubes=i, model=my_model)
            for t in i:
                for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
                    all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度

        # print("n = {0}, analysis time:{1}, {2}".format(self.n, time.time() - start, len(all_conc)))
        # all_conc = {}  # 记录结果
        # for t in tubes:
        #     for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
        #         all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度
        all_conc = sorted(all_conc.items(), key=lambda d: d[1], reverse=True)  # 排序

        error = {}  # 怀疑是错配的
        # k_cou = 0  # 记录浓度大于某个值的

        for k, v in all_conc:
            if k.count("+") == 2 and v > self.first_check:  # 将浓度换成输入的浓度
                # k_cou += 1
                # 根据：分割k， 然后根据名字具有顺序关系，然后确定是不是正确的配对
                tem_split = k.split('+')
                t1 = int(tem_split[0].split(':')[1][1:])
                t2 = int(tem_split[1].split(':')[1][1:])
                t3 = int(tem_split[2].split(':')[1][1:-1])
                set_t = {abs(t1 - t2), abs(t1 - t3), abs(t2 - t3)}

                if set_t != {1, self.len_g1, self.len_g1 - 1}:
                    tem_arr = [t1, t2, t3]
                    tem_arr.sort()
                    key = ''
                    for i in tem_arr:
                        if i >= self.len_g1:
                            key = key + 'R{0},'.format(i - self.len_g1)
                        else:
                            key = key + 'F{0},'.format(i)
                    key = key[:-1]
                    if key in error and v < error[key]:
                        continue
                    else:
                        error[key] = v

            elif v < self.first_check:  #
                break

        # second check
        # error_end = []  # 经过校验后还是错配的

        # 添加验证
        # print("验证:{0},{1},{2}".format(t1, t2, t3))
        # tem_err_list = self.get_tube(t1)
        # tem_err_list.extend(self.get_tube(t2))
        # tem_err_list.extend(self.get_tube(t3))
        # if self.verification_two(tem_err_list):
        #     error_end.append(k)

        # print("目标:{0},list1:{1},list2:{2}".format(int(self.len1 / 2), self.len_g1, self.len_g2))
        # print("出错的:{0},{1}".format(len(error), error))
        # print("最终检测还是出错的{0}".format(error_end))
        # print("浓度大于1e-9数:{0}".format(k_cou))
        # print('总数:{0}'.format(len(all_conc)))
        return error

    def get_strand_tube_all(self):
        # 获取试管中只有两条基因片段的所有情况

        self.list_g1 = self.list_g1[:-1]
        self.len_g1 = len(self.list_g1)
        strands = {}
        for i in range(0, self.len_g1):
            strands[Strand(self.list_g1[i], name="F{0}".format(i + 1))] = self.c_gene
            strands[Strand(self.list_g2[i], name="R{0}".format(i + 1))] = self.c_gene

        strands[Strand(self.F_Primer, name="F_Primer")] = self.c_gene
        strands[Strand(self.R_Primer, name="R_Primer")] = self.c_gene

        my_model = Model(material='dna', celsius=self.temp)
        t = Tube(strands=strands, complexes=SetSpec(max_size=2), name='t')  # complexes defaults to [A, B]
        tube_results = tube_analysis(tubes=[t], model=my_model)

        all_conc = {}
        for my_complex, conc in tube_results.tubes[t].complex_concentrations.items():
            all_conc[my_complex.name] = conc  # 反应后每个试管中DNA的浓度
        all_conc = sorted(all_conc.items(), key=lambda d: d[1], reverse=True)  # 排序

        # print(all_conc)

        error = {}
        for key, val in all_conc:
            if key.count('+') == 1 and val > self.first_check:
                tem_split = key[1:-1].split('+')  # # 先去除括号，在根据+分割字符串

                t1 = tem_split[0][1:]
                t2 = tem_split[1][1:]
                # print(t1, t2)
                # print(t1.isdigit(), t2.isdigit())
                if not t1.isdigit() or not t2.isdigit():
                    # 处理含有primer的片段，都是错配的
                    # primer+pirmer, primer+x         ,  x+primer
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val

                elif tem_split[0][0] == tem_split[1][0]:  # 错配
                    # Rx+Rx ; Fx+Fx
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                elif tem_split[0][0] == 'F' and int(t1) - int(t2) not in [0, 1]:
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
                elif tem_split[0][0] == 'R' and int(t1) - int(t2) not in [0, -1]:
                    error['{0}, {1}'.format(tem_split[0], tem_split[1])] = val
            # 判断条件有待改进
            elif val < self.first_check:  #
                break

        # print(error)

        return error

# if __name__ == '__main__':
#
#     data1 = ['TAAGCACCTGTAGGATCGTACAGGTTTACGCAAGAAAATGGTTTGTT', 'ATAGTCGAATAACACCGTGCGTGTTGACTATTTTACCTCTGGCGG', 'TGATATACTAGAGAAAGAGGAGAAATACTAGATGACCATGATTACGCCAAGCG', 'CGCAATTAACCCTCACTAAAGGGAACAAAAGCTGGAGCTCCACCG', 'CGGTGGCGGCAGCACTAGAGCTAGTGGATCCCCCGG', 'GCTGTAGAAATTCGATATCAAGCTTATCGATACCGTCGACCTCGAGGGG', 'GGGCCCGGTACCCAATTCGCCCTATAGTGAGTCGTATTACGCG', 'CGCTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACC', 'TGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTC', 'GCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCG', 'CCCTTCCCAACAGTTGCGCAGCCTGAATAATAACGCTGATAGTGCTA', 'TGTAGATCGCTACTAGAGCCAGGCATCAAATAAAACGAAAGGCTCAGTCG', 'AGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTC', 'CTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGT', 'TTATACAGTTCCTGAAGATAGATTAAGGCAC']
#
#     data2 = ['TGTACGATCCTACAGGTGCTTA', 'ACGCACGGTGTTATTCGACTATAACAAACCATTTTCTTGCGTAAACC', 'TCTAGTATTTCTCCTCTTTCTCTAGTATATCACCGCCAGAGGTAAAATAGTCAAC', 'TTCCCTTTAGTGAGGGTTAATTGCGCGCTTGGCGTAATCATGGTCA', 'AGTGCTGCCGCCACCGCGGTGGAGCTCCAGCTTTTG', 'CGATAAGCTTGATATCGAATTTCTACAGCCCGGGGGATCCACTAGCTC', 'CGAATTGGGTACCGGGCCCCCCCTCGAGGTCGACGGTA', 'ACGACGGCCAGTGAGCGCGCGTAATACGACTCACTATAGGG', 'GCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTG', 'CGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAG', 'GCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCT', 'CCTGGCTCTAGTAGCGATCTACACTAGCACTATCAGCGTTATTATTCAGGC', 'CAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGAT', 'AGGTGAGCCAGTGTGACTCTAGTAGAGAGCGTTCACCGACAAACAA', 'GTGCCTTAATCTATCTTCAGGAACTGTATAAACGCAGAAAGGCCCACCC']
#     data3 = 29
#     data4 = 37
#     data5 = 1e-8
#
#     veri = Verification(data1, data2[1:], data3, data4, data5)
#     ver_info = veri.analysis_two()
#     print(len(ver_info), ver_info)
#
#     ver_info1 = veri.analysis_three()
#     print(len(ver_info1), ver_info1)
