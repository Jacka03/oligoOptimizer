from django.contrib.staticfiles.storage import staticfiles_storage
from django.urls import reverse
from jinja2 import Environment


def environment(**options):
    env = Environment(**options)
    # env.globals.update({  # 修改
    #     'static': staticfiles_storage.url,  # 模板出现static，调用url函数
    #     'url': reverse  # 模板出现url，调用reverse函数
    # })
    return env
