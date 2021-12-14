import math

import numpy as np
from matplotlib import pyplot as plt


def show_w(x, y, head):
    plt.scatter(x, y)
    plt.plot(x, y)
    tem_str = "min:{:3f},max:{:3f},max-min={:3f},std={:4f}".format(min(y), max(y), max(y) - min(y), np.std(y))
    plt.legend([tem_str])
    plt.title(head)
    plt.show()


class Splicing:

    def __init__(self, input_info):
        self.input_info = input_info  # 各种离子信息
        self.gene = input_info['gene'].upper()  # 基因序列
        self.gene_len = input_info['geneLen']  # 基因序列
        self.res_type = input_info['resultType']  # 结果类型
        self.result = input_info['result']

        self.min_len = int(input_info['minLen'])
        self.max_len = int(input_info['maxLen']) + 1
        self.segment_len = self.max_len - self.min_len
        self.count = 20

        self.tail = False

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

        # c_Na = 1.3 # mmol /

        c_Na = self.input_info['Na']  #
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
        计算第三段到最后的切割位点
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
                    if sec_cut > self.gene_len - 1:
                        sec_cut = self.gene_len - 1
                    tem_tm = self.cal_tm(self.gene[fir_cut: sec_cut])  # 计算这段gene的tm
                    bef_tm = result[i, 1::2]  # 取出前面所有tm
                    bef_tm = np.append(bef_tm, tem_tm)  # 将这段gene的tm添加到之前中
                    if tm_mea != 0.:
                        bef_tm = np.append(bef_tm, tm_mea)  # 将这段gene的tm添加到之前中
                    tm_std = np.std(bef_tm)  # 计算标准差
                    bef_arr = result[i, :]  # 获取数组，转化为列表
                    bef_arr = bef_arr.tolist()
                    tem_gene_tm = [sec_cut, tem_tm, tm_std]
                    tem_list = bef_arr + tem_gene_tm

                    if fir_cut + self.min_len > self.gene_len - 1:
                        answer1_tem.append(tem_list)  # TODO最后一段是独立好还是分开好
                        break
                    elif sec_cut == self.gene_len - 1:
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

        show_w(index_res, tm_res, "greedy")
        return index_res, tm_res

    def cal(self):
        index, tm = self.cal_next_tm()
        index1, tm1 = self.cal_next_tm(np.mean(tm))


if __name__ == '__main__':
    data = {'email': '758168660@qq.com', 'geneLen': 638, 'result': 'res1', 'minLen': 20, 'maxLen': 30, 'resultType': 'Gap', 'verification': 'No', 'pools': 1, 'geneName': 'name', 'geneDesc': 'description', 'temperature': 37, 'concentrations': 1, 'gene': 'taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaataataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata', 'K': 50, 'Mg': 8, 'dNTPs': 4, 'Tris': 10, 'oligo': 10, 'primer': 400, 'Na': 1.2}
    sp = Splicing(data)
    sp.cal()