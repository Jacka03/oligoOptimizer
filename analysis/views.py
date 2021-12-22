import json

from django.http import JsonResponse, Http404
from django.shortcuts import render

# Create your views here.
from django.views import View

from analysis import models
from analysis.tool.splicing import Splicing
from analysis.tool.utils import is_dna, reverse_comple, preprocessing_data, get_res_info


class AssemblyView(View):

    def get(self, request):
        return render(request, 'assembly.html')

    def post(self, request):
        try:
            data = json.loads(request.body)
            data = preprocessing_data(data)

            # add to db
            # models.GeneInfo.objects.create(email=data['email'], gene_len=data['geneLen'],
            #                                pools=data['pools'], min_len=data['minLen'], max_len=data['maxLen'])

            splic = Splicing(data)
            next_cal, info = splic.cal()

            # add cal info to context
            tem_res = get_res_info(info)
            # return info
            context = {
                'info': info.get('result'),
                'resInfo': tem_res,
                'nextCal': next_cal
            }
            # print(context)
            # if data.get('verification') == 'Yes':
            #
            #     conc = data['concentrations'] * 1e-8
            #     # 分析过程
            #     analy = Analysis(next_cal[0], next_cal[1][1:], next_cal[2], data['temperature'], conc)
            #     analy_info = analy.analysis_two()
            #     analy_info.update(analy.analysis_three())
            #
            #     analy_info_list = []
            #     for key, value in analy_info.items():
            #         analy_info_list.append({
            #             'key': key,
            #             'value': value,
            #         })
            #     context['analyInfo'] = analy_info_list

            arr = [context]
            context = {'arr': arr}
        except Exception as err:
            # raise Http404('genenenenen')
            return JsonResponse({"error": err})

        return JsonResponse(context)


class AssemblyPoolsView(View):

    def post(self, request):
        # try:
        data = json.loads(request.body)
        data = preprocessing_data(data)

        pools = int(data['pools'])
        # print("pools:{0}".format(pools))
        splic = Splicing(data)
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

        context = {'arr': arr}
        # except Exception as e:
        #     print("error")
        #     context = {'error': 'error'}

        return JsonResponse(context)