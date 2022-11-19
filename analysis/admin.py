from django.contrib import admin

# Register your models here.
from . import models

admin.site.register(models.User)


class GeneInfoAdmin(admin.ModelAdmin):
    list_display = ('id', 'email', 'ip', 'ip_addr', 'date', 'gene_len', 'min_len', 'max_len', 'pools', 'get_avg_pool', 'assembly_time')

    @admin.display(ordering='get_avg_pool')
    def get_avg_pool(self, obj):
        return obj.gene_len / obj.pools

    list_per_page = 25


admin.site.register(models.GeneInfo, GeneInfoAdmin)


class VerificationAdmin(admin.ModelAdmin):
    list_display = ('id', 'cube_count', 'gene_segment_count',
                    'verification_two_time', 'verification_three_time', 'all_time')

    @admin.display(ordering='all_time')
    def all_time(self, obj):
        return obj.verification_two_time + obj.verification_three_time

    list_per_page = 25


admin.site.register(models.VerificationInfo, VerificationAdmin)
