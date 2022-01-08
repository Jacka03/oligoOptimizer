import math

import numpy as np
from matplotlib import pyplot as plt


def show_w(x, y, title):
    # 打印index与tm的图
    plt.scatter(x, y)
    plt.plot(x, y)
    tem_str = "min:{:3f},max:{:3f},max-min={:3f},std={:4f}".format(min(y), max(y), max(y) - min(y), np.std(y))
    plt.legend([tem_str])
    plt.title(title)
    plt.show()


def show_each_gene_len(index_list):
    # 打印每个片段的长度
    res = []
    for i in range(1, len(index_list)):
        res.append(index_list[i] - index_list[i - 1])
    print(res)
    print("min:{0}, max:{1}".format(min(res), max(res)))


class Splicing:

    def __init__(self, input_info):
        self.input_info = input_info  # 各种离子信息
        self.gene = input_info['gene']  # 基因序列
        self.gene_len_sor = len(self.gene)  # 基因序列
        self.gene_len = len(self.gene)  # 基因序列
        self.res_type = input_info['resultType']  # 结果类型
        self.result = input_info['result']

        self.min_len = int(input_info['minLen'])
        self.max_len = int(input_info['maxLen']) + 1
        self.segment_len = self.max_len - self.min_len
        self.count = 20

        self.tail = False
        self.above_tail = 0
        self.primer_min_len = 18
        self.primer_max_len = 22 + 1

    def cal_tm(self, temp_gene):
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

        H = AATT * (-7.6) + ATTA * (-7.2) + TAAT * (-7.2) + CAGT * (-8.5) + GTCA * (-8.4) + CTGA * (-7.8) + GACT * (
            -8.2) + CGGC * (-10.6) + GCCG * (-9.8) + GGCC * (-8.0) + 0.2 + 2.2
        S = AATT * (-21.3) + ATTA * (-20.4) + TAAT * (-21.3) + CAGT * (-22.7) + GTCA * (-22.4) + CTGA * (
            -21.0) + GACT * (
                -22.2) + CGGC * (-27.2) + GCCG * (-24.4) + GGCC * (-19.9) - 5.7 + 6.9 - 1.4

        # TODO 钠离子浓度是多少？

        # c_Na = self.input_info['Na']  #
        c_Na = 1.3  # mmol /
        #
        c_K = self.input_info['K'] / 1000  #
        c_Mg = self.input_info['Mg'] / 1000  #
        c_dNTPs = self.input_info['dNTPs'] / 1000
        c_Tris = self.input_info['Tris'] / 1000  # mol / L#

        c_oligo = self.input_info['oligo'] / 1e9  # 寡核苷酸
        c_t = self.input_info['primer'] / 1e9  # 引物cd

        # TODO 当premer不是远大于oligo时，c_t需要重新计算

        # TODO Mon离子浓度初始化
        # c_Mon = 0.005
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

    def cal_first_tm(self, mean_tm=0.):
        """
        计算第一段与第二段的tm值并且计算其标准差
        :param mean_tm:
        :return:
        """
        result = []
        for i in range(self.segment_len):
            mid_cut = self.min_len + i
            fir_tm = self.cal_tm(self.gene[:mid_cut])
            if mean_tm == 0.:
                for j in range(self.segment_len):
                    end_cut = mid_cut + self.min_len + j
                    sec_tm = self.cal_tm(self.gene[mid_cut:end_cut])
                    result.append([mid_cut, fir_tm, end_cut, sec_tm, np.std([fir_tm, sec_tm])])
            else:  # 使用经验值作为开始
                result.append([mid_cut, fir_tm, np.std([mean_tm, fir_tm])])
        return result

    def choose(self, tm_list, count=1):
        """
        根据tm标准差降序排序，返回前count个标准差小的分割方法
        :param tm_list:
        :param count:
        :return:
        """
        tm_list = np.array(tm_list)
        tm_list = tm_list[np.lexsort(tm_list.T)]
        # count = 10  # 在二维数组中获取前count个tm标准差最小的
        if len(tm_list) > count:
            return tm_list[:count]
        return tm_list

    def cal_next_tm(self, tm_mea=0.):
        """
        计算第三段到最后的切割位点,贪心算法
        :return:
        """
        result = self.cal_first_tm(tm_mea)
        result = self.choose(result, self.count)
        result = np.delete(result, -1, axis=1)  # 删除最后一列

        fir_ans_tem = []  # 初步切割结果
        answer1_tem = []
        # 尝试控制answer不为0
        while len(fir_ans_tem) == 0:
            tem_res = []
            for i in range(len(result)):  # 遍历上一轮选择到的最优的
                fir_cut = int(result[i, -2])  # 这段gene开始
                # TODO 如果这个长度小鱼最小切割位点，退出（推出前应该如何）
                if fir_cut + self.min_len > self.gene_len:
                    fir_ans_tem = result
                    break
                for j in range(self.segment_len):  #
                    sec_cut = fir_cut + self.min_len + j  # 这段gene结束
                    if sec_cut > self.gene_len:
                        sec_cut = self.gene_len
                    tem_tm = self.cal_tm(self.gene[fir_cut: sec_cut])  # 计算这段gene的tm
                    bef_tm = result[i, 1::2]  # 取出前面所有tm
                    bef_tm = np.append(bef_tm, tem_tm)  # 将这段gene的tm添加到之前中
                    # if tm_mea != 0.:
                    #     bef_tm = np.append(bef_tm, tm_mea)  # 将这段gene的tm添加到之前中
                    tm_std = np.std(bef_tm)  # 计算标准差
                    bef_arr = result[i, :]  # 获取数组，转化为列表
                    bef_arr = bef_arr.tolist()
                    tem_gene_tm = [sec_cut, tem_tm, tm_std]
                    tem_list = bef_arr + tem_gene_tm

                    if fir_cut + self.min_len > self.gene_len:
                        answer1_tem.append(tem_list)  # TODO最后一段是独立好还是分开好
                        break
                    elif sec_cut == self.gene_len:
                        fir_ans_tem.append(tem_list)
                        break
                    else:
                        tem_res.append(tem_list)
            # 可能刚刚好处理完
            if len(tem_res) != 0:
                tem_res = self.choose(tem_res, self.count)
                result = np.delete(tem_res, -1, axis=1)  # 删除最后一列

        # 挑选结果
        fir_ans_tem = self.choose(fir_ans_tem)
        if len(answer1_tem) > 0 and len(answer1_tem[0]) == len(answer1_tem[-1]):
            answer1_tem = self.choose(answer1_tem)
            if fir_ans_tem[0, -1] > answer1_tem[0, -1]:
                fir_ans_tem = answer1_tem
        index_res = np.array(fir_ans_tem[0][:-1:2])
        tm_res = np.array(fir_ans_tem[0][1::2])

        # show_w(index_res, tm_res, "greedy")
        return index_res, tm_res

    def cal_all_tm(self, arr):
        """
        求这种切割位点的tm的标准差
        :param arr: 整个基因片段的切割位点
        :return: np.std(tm_list)：标准差， tm_list：每段的tm组成的list
        """
        tm_list = []
        arr = arr.astype(int)
        for i in range(1, len(arr)):
            tm_t = self.cal_tm(self.gene[arr[i - 1]: arr[i]])
            tm_list.append(tm_t)
        return np.std(tm_list), tm_list

    def iteration(self, index_list, tm_list):  # 全局迭代，从左到右
        """
        全局迭代，
        :param index_list:上次得到最好的剪切位置
        :param tm_list:每个剪切片段的tm
        :return:
        """
        flag = [-1, -2, -3]  # 提前停止遍历的终止判断
        tem_max_len = self.max_len + 5
        tem_min_len = self.min_len - 5
        bias_len = 5
        while flag[-1] != flag[-2]:
            for i in range(len(tm_list)):  # 对于整条
                left = index_list[i]
                right = index_list[i + 1]
                # 遍历
                tem_result = []
                for j in range(1 - bias_len, bias_len):  # 对于每个片段,
                    tem_left = left + j  # 左边新的切割位点
                    if left == 0 and tem_left != 0:
                        # 第一个切割位点不能动
                        continue

                    # if tem_left < 0:
                    #     # 当第一个位点小于0时，舍弃该情况
                    #     continue
                    # 左边移动后前一个片段的长度判断：
                    if i >= 1 and (
                            tem_left - index_list[i - 1] < tem_min_len or tem_left - index_list[i - 1] > tem_max_len):
                        # 当左边位点移动导致左边片段过长或者过短时，舍弃该情况
                        continue

                    for k in range(1 - bias_len, bias_len):
                        tem_right = right + k  # 右边新的切割位点

                        if right == self.gene_len and tem_right != self.gene_len:
                            # 原始序列最后一个
                            continue

                        # if tem_right >= self.gene_len:
                        #     # 右边切割位点超越右边界
                        #     continue
                        # 右节点移动后下一个片段长度判断
                        if i + 1 < len(tm_list) and (
                                index_list[i + 2] - tem_right > tem_max_len
                                or index_list[i + 2] - tem_right < tem_min_len):
                            # 当右边位点移动导致右边片段过长或者过短时，舍弃该情况
                            continue
                        # 当前片段长度是否正确
                        if tem_min_len > (tem_right - tem_left) or (tem_right - tem_left) > tem_max_len:
                            # 当前片段过长或者过短时
                            continue

                        tem_index_list = index_list.copy()
                        tem_index_list[i] = tem_left
                        tem_index_list[i + 1] = tem_right
                        tem_std, _ = self.cal_all_tm(tem_index_list)
                        tem_result.append([tem_left, tem_right, tem_std])

                tem_result = self.choose(tem_result, count=1)
                index_list[i] = tem_result[0, 0]
                index_list[i + 1] = tem_result[0, 1]
                best_std, tm_list = self.cal_all_tm(index_list)  # 本次迭代得到最好的std和tm_list

            flag.append(best_std)
        #     show_w(index_list[1:], tm_list, "d")
        # show_w(index_list[1:], tm_list, "iteration")

        return index_list, tm_list

    def overlap(self, index_list, tm_list):
        # 动态规划
        index_list = [int(i) for i in index_list]
        gene_list = []
        for i in range(len(tm_list)):  # 将gene截取出来存放在一个二维list中
            # 【原来第一个切割位点，修改后第一个切割位点，原来第二个切割位点，修改后第二个切割位点，修改后片段tm】
            gene_list.append([index_list[i], index_list[i], index_list[i + 1], index_list[i + 1], tm_list[i]])

        over_size = 6
        temp_avg_tm = [0, 1]  # flag 结束迭代的标志

        tem_max_op = 10  # 相邻两个片段间隔最大值

        tem_min_len = self.min_len - 10  # 切割后每个片段最小值

        x = 0
        while temp_avg_tm[-1] != temp_avg_tm[-2]:  # 迭代，终止条件

            for i in range(len(gene_list)):
                new_tm_list = tm_list.copy()
                tem_result = []
                for j in range(over_size):
                    left = gene_list[i][1] + j
                    if gene_list[i][1] == 0 and left != 0:
                        # 第一个位点不处理
                        continue

                    for k in range(over_size):

                        right = gene_list[i][2] - k
                        if gene_list[i][2] == self.gene_len and right != self.gene_len:
                            continue

                        if right - left < tem_min_len:  # 长度小于限定值
                            continue
                        if (i == 0 and gene_list[i][1] - 0 > tem_max_op) or (
                                i > 0 and left - gene_list[i - 1][2] > tem_max_op):  # 这段与前一段的距离
                            continue
                        if (i + 1 == len(gene_list) and len(self.gene) - gene_list[i][2] - 1 > tem_max_op) or (
                                i + 1 < len(gene_list) and gene_list[i + 1][1] - right > tem_max_op):
                            continue

                        tem_tm = self.cal_tm(self.gene[left: right])
                        new_tm_list[i] = tem_tm
                        tm_std = np.std(new_tm_list)
                        tem_result.append([left, right, tem_tm, tm_std])
                if len(tem_result) == 0:
                    continue
                tem_result = self.choose(tem_result, 1)
                gene_list[i][1] = int(tem_result[0, 0])
                gene_list[i][2] = int(tem_result[0, 1])
                tm_list[i] = tem_result[0, 2]
                temp_avg_tm.append(tem_result[0, 3])
            x = x + 1
            # show_w(index_list[1:], tm_list, x)

        # show_w(index_list[1:], tm_list, "overlap")

        a = np.argsort(tm_list)

        # 经过剪切后，在迭代一次，进行扩展
        for i in range(len(a)):
            if a[i] < 1 or a[i] + 2 > len(a):  # 先不管第一段和最后一段
                continue
            test_tm_list = tm_list.copy()
            test_result = []
            for j in range(int(gene_list[a[i]][1] - gene_list[a[i] - 1][2])):
                for k in range(int(gene_list[a[i] + 1][1] - gene_list[a[i]][2])):
                    test_gene = self.gene[gene_list[a[i]][1] - j: gene_list[a[i]][2] + k]
                    test_tm = self.cal_tm(test_gene)
                    test_tm_list[a[i]] = test_tm
                    test_std = np.std(test_tm_list)
                    test_result.append([gene_list[a[i]][1] - j, gene_list[a[i]][2] + k, test_tm, test_std])
            if len(test_result) == 0:
                continue
            # 下表从0开始
            test_result = self.choose(test_result, 1)
            gene_list[a[i]][1] = test_result[0, 0]
            gene_list[a[i]][2] = test_result[0, 1]
            gene_list[a[i]][4] = test_result[0, 2]
            tm_list[a[i]] = test_result[0, 2]

        # show_w(index_list[1:], tm_list, "end")
        return gene_list

    def input_tail(self, tem_index, tem_tm):
        if not isinstance(tem_index, list):
            tem_index = tem_index.tolist()
        if not isinstance(tem_tm, list):
            tem_tm = tem_tm.tolist()

        str = "CAGTTCCTGAAGATAGATTAAGGCACCGTGATGAACGTATGCACAGCTTCCGAAGGTGAGCCAGTGTGACTCTAAACAGATAAAACGAAAG"
        tem_gene = self.gene + str

        tem_ans = []
        for i in range(self.segment_len):
            end_cut = int(tem_index[-1] + self.min_len + i)
            end_tm = self.cal_tm(tem_gene[int(tem_index[-1]): end_cut])
            # cal
            tm_list = tem_tm + [end_tm]
            end_std = np.std(tm_list)
            tem_ans.append([end_cut, end_tm, end_std])

        # print(len(self.gene), self.gene_len, self.gene_len_sor, len(tem_gene))
        # print(tem_ans)
        tem_res = self.choose(tem_ans, 1)
        tem_index.append(tem_res[0][0])
        tem_tm.append(tem_res[0][1])

        self.gene = tem_gene
        self.gene_len = len(self.gene)
        self.tail = True
        return tem_index, tem_tm

    def cal_mean_std(self, index):
        # 根据一个切割位点list转化为[[i1, i1, i2, i2, tm], [i2, i2, i3, i3, tm], ...]的形式, 方便后面计算这种切割的overlap的tm
        tem_res = []
        for i in range(1, len(index)):
            tem_gene = self.gene[int(index[i - 1]):int(index[i])]
            tem_res.append([int(index[i - 1]), int(index[i - 1]), int(index[i]), int(index[i]), self.cal_tm(tem_gene)])
        return tem_res

    def return_result(self, index, tm):
        if self.res_type == 'Gapless':
            return self.cal_mean_std(index)
        elif self.res_type == 'Gap':
            return self.overlap(index, tm)

    def id_dna_complete(self, cut_of_index):
        # 查看dna是否是完整的
        index = int(cut_of_index[-1][2])
        while index < self.gene_len_sor:
            self.above_tail = self.above_tail + 1
            tem_tm = []
            for i in cut_of_index:
                tem_tm.append(i[-1])
            # if index != self.gene_len:
            self.gene = self.gene + "CAGTTCCTGAAGATAGATTAAGGCACCGTGATGAACGTATGCACAG"
            self.gene_len = len(self.gene)
            self.tail = True

            tem_ans = []
            for i in range(self.segment_len):
                end_cut = index + self.min_len + i
                end_tm = self.cal_tm(self.gene[index: end_cut])
                # cal
                tm_list = tem_tm + [end_tm]
                end_std = np.std(tm_list)
                tem_ans.append([end_cut, end_tm, end_std])

            tem_res = self.choose(tem_ans, 1)
            cut_of_index.append([index, index, tem_res[0][0], tem_res[0][0], tem_res[0][1]])
            index = int(cut_of_index[-1][2])

        return cut_of_index

    def add_tail(self, cut_of_index):
        self.above_tail = self.above_tail + 1

        index = int(cut_of_index[-1][2])
        tem_tm = []
        for i in cut_of_index:
            tem_tm.append(i[-1])

        str = "CAGTTCCTGAAGATAGATTAAGGCACCGTGATGAACGTATGCACAGCTTCCGAAGG"
        self.gene = self.gene + str
        self.gene_len = len(self.gene)
        self.tail = True

        tem_ans = []
        for i in range(self.segment_len):
            end_cut = index + self.min_len + i
            end_tm = self.cal_tm(self.gene[index: end_cut])
            # cal
            tm_list = tem_tm + [end_tm]
            end_std = np.std(tm_list)
            tem_ans.append([end_cut, end_tm, end_std])

        tem_res = self.choose(tem_ans, 1)
        cut_of_index.append([index, index, tem_res[0][0], tem_res[0][0], tem_res[0][1]])

        return cut_of_index

    def get_primer(self, oligo, data):
        # print(data)
        result = []
        for i in range(self.primer_min_len, self.primer_max_len):
            # 当oligo长度过长时
            tem_data = data + [self.cal_tm(oligo[:i])]
            result.append([i, np.std(tem_data)])
        # print(result)
        result = self.choose(result)
        oligo = oligo[:int(result[0][0])]

        return oligo

    def get_oligo(self, index_list):
        # print(len(index_list))
        # 将两个overlap拼接成一个Oligo
        dnaTable = {
            "A": "T", "T": "A", "C": "G", "G": "C"
        }
        gene_complement = ""  # DNA互补链
        for ele in self.gene:
            gene_complement += dnaTable[ele]

        # 存储DNA双链上面的oligo

        oligo_list = []
        above_oligo_info = []
        for i in range(0, len(index_list) - 1, 2):
            oligo = self.gene[int(index_list[i][1]): int(index_list[i + 1][2])]
            index = int(i / 2)
            oligo_list.append(['F{0}'.format(index), oligo])
            overlap = self.gene[int(index_list[i][1]): int(index_list[i][2])]
            above_oligo_info.append(
                ['F{0}'.format(index), oligo, round(self.cal_tm(overlap), 2), len(overlap), len(oligo)])

        oligo = self.gene[int(index_list[0][1]): int(index_list[0][2])]
        oligo_list.append(['F_Primer', oligo])  # 用于验证
        above_oligo_info.append(['F_Primer', oligo, round(self.cal_tm(oligo), 2), '', len(oligo)])

        # 存储DNA双链下面的oligo
        below_oligo_info = []
        for i in range(1, len(index_list) - 1, 2):
            oligo = gene_complement[int(index_list[i][1]): int(index_list[i + 1][2])][::-1]  # 转化成从左到右是5'到3'
            index = int((i - 1) / 2)
            oligo_list.append(['R{0}'.format(index), oligo])
            overlap = gene_complement[int(index_list[i][1]): int(index_list[i][2])][::-1]
            below_oligo_info.append(
                ['R{0}'.format(index), oligo, round(self.cal_tm(overlap), 2), len(overlap), len(oligo)])
        # 最后一片
        oligo = gene_complement[int(index_list[-1][1]): int(index_list[-1][2])][::-1]
        if self.tail:  # R_Primer避开添加上去的tail
            oligo = gene_complement[self.gene_len_sor - self.max_len: self.gene_len_sor][::-1]
        oligo_list.append(['R_Primer', oligo])
        below_oligo_info.append(['R_Primer', oligo, round(self.cal_tm(oligo), 2), '', len(oligo)])

        result = []  # [[Label, oligo, tm, overlap_len, oligo_len], ]
        for i in range(len(above_oligo_info)):
            result.append(above_oligo_info[i])
            result.append(below_oligo_info[i])

        data = []
        for i in range(len(result) - 2):
            data.append(result[i][2])

        # 只有得到data数组才能继续
        if result[-2][-1] >= self.primer_max_len:
            oligo = self.get_primer(result[-2][1], data)
            oligo_list[len(above_oligo_info)-1] = ['F_Primer', oligo]
            result[-2] = ['F_Primer', oligo, round(self.cal_tm(oligo), 2), '', len(oligo)]

        if result[-1][-1] >= self.primer_max_len:
            oligo = self.get_primer(result[-1][1], data)
            oligo_list[-1] = ['R_Primer', oligo]
            result[-1] = ['R_Primer', oligo, round(self.cal_tm(oligo), 2), '', len(oligo)]

        # overlap的信息
        overlap_data = {
            'min': round(min(data), 2),
            'max': round(max(data), 2),
            'range': round(max(data) - min(data), 2),
            'mean': round(np.mean(data), 2),
            'std': round(np.std(data), 2),

            'result': result,
        }
        # self.tail = False
        if self.above_tail > 2:
            print(self.gene, self.min_len, self.max_len)
            raise Exception

        if self.above_tail == 2:
            print(self.gene[self.gene_len_sor: int(index_list[-2][2])])

        overlap_data['tail_reverse'] = ''
        if self.tail:
            # print(self.gene_len_sor, cut_of_index[-1][2])
            # print(index_list[-1][2])
            overlap_data['tail'] = self.gene[self.gene_len_sor: int(index_list[-1][2])]
            overlap_data['tail_reverse'] = gene_complement[self.gene_len_sor: int(index_list[-1][2])][::-1]

        return overlap_data, oligo_list

    def cal(self):
        index, tm = self.cal_next_tm()
        # show_w(index, tm, "1")
        # index, tm = self.cal_next_tm(np.mean(tm))

        index = np.insert(index, 0, [0])
        index, tm = self.iteration(index, tm)

        # Gap or Gapless
        cut_of_index = self.return_result(index, tm)

        # 查看是否包含完整的链进去了
        cut_of_index = self.id_dna_complete(cut_of_index)
        # add tail
        if len(cut_of_index) % 2 == 0:
            cut_of_index = self.add_tail(cut_of_index)

        info, oligo = self.get_oligo(cut_of_index)
        # info for valification
        return oligo, info

    def cal_for_pool(self):
        _, tm = self.cal_next_tm()
        index, tm = self.cal_next_tm(np.mean(tm))
        index = [int(i) for i in index]
        return index, tm  # TODO

    def cal_for_each_pool(self, index, tm):
        # add tail
        if len(index) % 2 == 0:
            # print(len(index), index[-1])
            index, tm = self.input_tail(index, tm)
            # print(len(index), index[-1])

        # 对整体遍历
        index = np.insert(index, 0, [0])
        index, tm = self.iteration(index, tm)

        cut_of_index = self.return_result(index, tm)

        info, oligo = self.get_oligo(cut_of_index)
        # info for valification
        return oligo, info


if __name__ == '__main__':
    data = {'email': '758168660@qq.com', 'geneLen': 638, 'result': 'res1', 'minLen': 20, 'maxLen': 30,
            'resultType': 'Gap', 'verification': 'No', 'pools': 1, 'geneName': 'name', 'geneDesc': 'description',
            'temperature': 37, 'concentrations': 1, 'gene': 'atgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaataataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctt', 'K': 50, 'Mg': 8, 'dNTPs': 4, 'Tris': 10, 'oligo': 10, 'primer': 400, 'Na': 1.2}
    data['gene'] = data['gene'].upper()
    sp = Splicing(data)
    next_cal, info = sp.cal()
    print(info)
    # taagc