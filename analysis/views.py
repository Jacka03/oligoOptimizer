import io
import json
import re
import time

import pandas as pd

from django.http import JsonResponse, Http404, HttpResponse
from django.shortcuts import render
from django.views import View

from analysis import models
from analysis.tool.splicing import Splicing
from analysis.tool.utils import preprocessing_data, get_res_info, get_ip
from analysis.tool.verification import Verification


class AssemblyView(View):

    def get(self, request):
        return render(request, 'assembly.html')

    def post(self, request):
        try:
            user_ip = request.META.get('REMOTE_ADDR')
            data = json.loads(request.body)
            data = preprocessing_data(data)
            ip_data = get_ip(user_ip)

            splic = Splicing(data)
            start = time.time()
            next_cal, info = splic.cal()
            end = time.time()
            # add to db
            models.GeneInfo.objects.create(email=data['email'], gene_len=data['geneLen'], pools=data['pools'],
                                           ip=ip_data['ip'], ip_addr=ip_data['address'],
                                           min_len=data['minLen'], max_len=data['maxLen'], assembly_time=end - start)

            models.Data.objects.create(email=data['email'], gene=data['gene'],
                                       min_len=data['minLen'], max_len=data['maxLen'], pools=data['pools'],
                                       k=data['K'], mg=data['Mg'], dntps=data['dNTPs'], tris=data['Tris'],
                                       primer=data['primer'], oligo=data['oligo'],
                                       gap=data['resultType'], verification=data['verification'])

            # add cal info to context
            tem_res = get_res_info(info)
            # return info
            context = {
                'info': info.get('result'),
                'resInfo': tem_res,
                'nextCal': next_cal,
                'tail_reverse': info.get('tail_reverse'),
                'above_tail': info.get('above_tail')
            }
            # print(context)
            if data.get('verification') == 'Yes':
                oligo_conc = data['oligoConc'] * 1e-9
                primer_conc = data['primerConc'] * 1e-9

                if primer_conc == 0 or primer_conc == 0.:
                    print(primer_conc)
                    primer_conc = oligo_conc

                # 分析过程
                analy = Verification(next_cal, data['temperature'], oligo_conc, primer_conc)

                analy_info = analy.get_strand_tube_all()

                analy_info_list = []
                for key, value in analy_info.items():
                    # analy_info_list.append({
                    #     'key': key,
                    #     'value': value,
                    # })
                    analy_info_list.append([key, value[0], value[1]])
                context['analyInfo'] = analy_info_list

            arr = [context]
            context = {'arr': arr}
        except Exception as err:
            print(err)
            raise Http404(err)

        return JsonResponse(context)


class AssemblyPoolsView(View):

    def post(self, request):
        try:
            user_ip = request.META.get('REMOTE_ADDR')
            ip_data = get_ip(user_ip)

            data = json.loads(request.body)
            data = preprocessing_data(data)

            pools = int(data['pools'])
            # print("pools:{0}".format(pools))
            splic = Splicing(data)
            start = time.time()
            index, tm = splic.cal_for_pool()

            # overlap of each pool
            each_pool = round((len(index) - 1) / pools)
            each_pool = [each_pool + 1, each_pool + 1] if each_pool % 2 == 0 else [each_pool, each_pool + 2]

            # each_pool = each_pool + 1 if each_pool % 2 == 0 else each_pool
            # print("after pools:{0}".format(each_pool))

            gene = data['gene']
            arr = []
            left = 0
            for i in range(pools):
                right = left + each_pool[i % 2]

                if i == 0:
                    index_list = index[left: right]
                    gene_list = gene[0: index_list[-1]]
                    tm_list = tm[left: right]

                elif i == pools - 1:
                    index_list = index[left: right]
                    gene_list = gene[index[left - 1]: len(gene)]
                    index_list = [x - index[left - 1] for x in index_list]
                    tm_list = tm[left: right]
                elif i > 0:
                    index_list = index[left: right]
                    gene_list = gene[index[left - 1]: index_list[-1]]
                    index_list = [x - index[left - 1] for x in index_list]
                    tm_list = tm[left: right]

                left = right - 1

                data['gene'] = gene_list

                splic = Splicing(data)
                # print(len(index_list))
                next_cal, info = splic.cal_for_each_pool(index_list, tm_list)

                tem_res = get_res_info(info)

                context = {
                    'info': info.get('result'),
                    'resInfo': tem_res,
                    'nextCal': next_cal,
                    'tail_reverse': info.get('tail_reverse'),
                    'above_tail': info.get('above_tail'),

                    'temperature': data['temperature'],
                    'oligoConc': data['oligoConc'],
                    'primerConc': data['primerConc']
                }
                arr.append(context)
            end = time.time()
            models.GeneInfo.objects.create(email=data['email'], gene_len=data['geneLen'], pools=data['pools'],
                                           ip=ip_data['ip'], ip_addr=ip_data['address'],
                                           min_len=data['minLen'], max_len=data['maxLen'], assembly_time=end - start)

            models.Data.objects.create(email=data['email'], gene=data['gene'],
                                       min_len=data['minLen'], max_len=data['maxLen'], pools=data['pools'],
                                       k=data['K'], mg=data['Mg'], dntps=data['dNTPs'], tris=data['Tris'],
                                       primer=data['primer'], oligo=data['oligo'],
                                       gap=data['resultType'], verification=data['verification'])

            context = {'arr': arr}
        except Exception as err:
            print(err)
            raise Http404(err)

        return JsonResponse(context)


class AnalysisView(View):

    def get(self, request):
        return render(request, 'assembly.html')

    def post(self, request):
        try:
            next_cal = json.loads(request.body)
            # next_cal = data['nextCal']
            # 分析过程
            oligo_conc = next_cal[2] * 1e-9
            primer_conc = next_cal[3] * 1e-9

            if primer_conc == 0 or primer_conc == 0.:
                print(primer_conc)
                primer_conc = oligo_conc

            analy = Verification(next_cal[0], next_cal[1], oligo_conc, primer_conc)
            analy_info = analy.get_strand_tube_all()

            analy_info_list = []
            context = {}
            for key, value in analy_info.items():
                # analy_info_list.append({
                #     'key': key,
                #     'value': value,
                # })
                analy_info_list.append([key, value[0], value[1]])

            context['analyInfo'] = analy_info_list

        except Exception as err:
            print(err)
            raise Http404(err)

        return JsonResponse(context)


class DownloadView(View):  # 导出excel数据
    def get(self, request):
        # print("test success")
        return HttpResponse("get")

    def replace_table(self, tem):
        if len(tem[-1]['resInfo']) == 6:
            # print(tem[-1]['info'][-3][1])
            # print(re.sub(r'<(.*?)>', r'', tem[-1]['info'][-3][1], count=2))
            # tem = re.sub(r'<(.*?)>', r'', tem[-1]['info'][-3][1], count=2)
            tem[-1]['info'][-3][1] = re.sub(r'<(.*?)>', r'', tem[-1]['info'][-3][1], count=2)
            tem[-1]['info'][-4][1] = re.sub(r'<(.*?)>', r'', tem[-1]['info'][-4][1], count=2)
            # tem[-1]['info'][-3][1] == tem[-1]['info'][-3][1].replace()
        return tem

    # 将dataform转换成django-excel下载是的sheet
    def post(self, request):
        try:
            tem = json.loads(request.body)
            tem = self.replace_table(tem)
            output = io.BytesIO()  # 配置一个BytesIO 这个是为了转二进制流
            count = 0
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                for i in range(len(tem)):
                    # pool
                    data = pd.DataFrame(data=['pool:{0}'.format(i + 1)])
                    data.to_excel(writer, startrow=count, index=False, header=None)
                    count += 1

                    # analysis result
                    res_info = tem[i]['resInfo']
                    data = pd.DataFrame(data=res_info)
                    data.to_excel(writer, startrow=count, startcol=6, index=False)

                    if 'analyInfo' in tem[i]:
                        # analysis result
                        analy_info = tem[i]['analyInfo']
                        data = pd.DataFrame(data=analy_info)
                        data.to_excel(writer, startrow=count, startcol=9, index=False)

                    # result
                    info = tem[i]['info']
                    data = pd.DataFrame(data=info, columns=['index', 'gene', 'tm', 'overlap', 'length'])
                    data.to_excel(writer, startrow=count, index=False)
                    count += data.shape[0] + 2

            # data.to_excel(output, index=False)  # index=False 是为了不建立索引
            # append_df_to_excel(data, output, sheet_name='Sheet1', startcol=10,startrow=0,index=False)
            output.seek(0)  # 把游标归0

            response = HttpResponse()  # 创建一个HttpResponse
            # response["Content-Type"] = "application/vnd.ms-excel"  # 类型
            file_name = 'comment.xlsx'  # 文件名称 自定义
            response['Content-Type'] = 'application/octet-stream'
            response['Content-Disposition'] = 'attachment;filename="{0}"'.format(file_name)
            response.write(output.getvalue())  # 写入数据
            output.close()  # 关闭
        except Exception as err:
            print(err)
            raise Http404(err)
        return response  # 返回


class IntroductionView(View):
    def get(self, request):
        return render(request, 'introduction.html')

    def post(self, request):
        return


class DocView(View):
    def get(self, request):
        return render(request, 'doc.html')

    def post(self, request):
        return
