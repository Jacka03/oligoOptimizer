import json
import time
import io
import re

import pandas as pd

from django.http import JsonResponse, Http404, HttpResponse
from django.shortcuts import render, redirect

# Create your views here.
from django.views import View

from analysis import models
from analysis.tool.splicing import Splicing
from analysis.tool.utils import is_dna, reverse_comple, preprocessing_data, get_res_info
from analysis.tool.verification import Verification


class AssemblyView(View):

    def get(self, request):
        return render(request, 'assembly.html')

    def post(self, request):
        try:
            data = json.loads(request.body)
            data = preprocessing_data(data)

            start = time.time()
            splic = Splicing(data)
            next_cal, info = splic.cal()
            end = time.time()
            # add to db
            models.GeneInfo.objects.create(email=data['email'], gene_len=data['geneLen'], pools=data['pools'],
                                           min_len=data['minLen'], max_len=data['maxLen'], assembly_time=end-start)

            # add cal info to context
            tem_res = get_res_info(info)
            # return info
            context = {
                'info': info.get('result'),
                'resInfo': tem_res,
                'nextCal': next_cal
            }
            # print(context)
            if data.get('verification') == 'Yes':
                conc = data['concentrations'] * 1e-8
                # 分析过程
                analy = Verification(next_cal[0], next_cal[1][1:], next_cal[2], data['temperature'], conc, next_cal[3])

                start = time.time()
                analy_info = analy.get_strand_tube_all()

                end = time.time()
                # print(end - start)
                # models.VerificationInfo.objects.create(cube_count=8000, gene_segment_count=next_cal[2],
                #                                        verification_two_time=mid-start, verification_three_time=end-mid)
                # print("two:{0}, three:{1}".format(mid-start, end-mid))

                analy_info_list = []
                for key, value in analy_info.items():
                    analy_info_list.append({
                        'key': key,
                        'value': value,
                    })
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
            data = json.loads(request.body)
            data = preprocessing_data(data)

            pools = int(data['pools'])
            # print("pools:{0}".format(pools))
            splic = Splicing(data)
            start = time.time()

            index, tm = splic.cal_for_pool()

            # overlap of each pool
            each_pool = int(len(index) / pools)
            each_pool = each_pool + 1 if each_pool % 2 == 0 else each_pool
            # print("after pools:{0}".format(each_pool))

            gene = data['gene']
            # print("gene:{0}".format(len(gene)))
            arr = []
            for i in range(pools):
                if i == 0:
                    index_list = index[i * each_pool: (i + 1) * each_pool]
                    gene_list = gene[0: index_list[-1]]
                    tm_list = tm[i * each_pool: (i + 1) * each_pool]
                elif i == pools - 1:
                    index_list = index[(pools - 1) * each_pool: -1]
                    gene_list = gene[index[(pools - 1) * each_pool - 1]: len(gene)]
                    index_list = [x - index[(pools - 1) * each_pool - 1] for x in index_list]
                    tm_list = tm[(pools - 1) * each_pool: -1]
                elif i > 0:
                    index_list = index[i * each_pool: (i + 1) * each_pool]
                    gene_list = gene[index[i * each_pool - 1]: index_list[-1]]
                    index_list = [x - index[i * each_pool - 1] for x in index_list]
                    tm_list = tm[i * each_pool: (i + 1) * each_pool]

                data['gene'] = gene_list

                splic = Splicing(data)
                # print(len(index_list))
                next_cal, info = splic.cal_for_each_pool(index_list, tm_list)

                tem_res = get_res_info(info)

                context = {
                    'info': info.get('result'),
                    'resInfo': tem_res,
                    'nextCal': next_cal,
                    'temperature': data['temperature'],
                    'concentrations': data['concentrations']
                }
                arr.append(context)
            end = time.time()
            models.GeneInfo.objects.create(email=data['email'], gene_len=data['geneLen'], pools=data['pools'],
                                           min_len=data['minLen'], max_len=data['maxLen'], assembly_time=end-start)

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
            temp = next_cal[5] * 1e-8
            # print(temp)
            analy = Verification(next_cal[0], next_cal[1][1:], next_cal[2], next_cal[4], temp, next_cal[3])

            start = time.time()

            analy_info = analy.get_strand_tube_all()

            # analy_info = analy.analysis_two()
            # mid = time.time()
            # analy_info.update(analy.analysis_three())
            end = time.time()
            # print(end - start)

            # models.VerificationInfo.objects.create(cube_count=8000, gene_segment_count=next_cal[2],
            #                                        verification_two_time=mid-start, verification_three_time=end-mid)

            analy_info_list = []
            context = {}
            for key, value in analy_info.items():
                analy_info_list.append({
                    'key': key,
                    'value': value,
                })
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