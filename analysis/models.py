from django.db import models


# Create your models here.

class GeneInfo(models.Model):
    email = models.EmailField(verbose_name='email')
    date = models.DateTimeField(auto_now_add=True)
    gene_len = models.IntegerField(verbose_name='gene_len')
    pools = models.IntegerField(verbose_name='pools')
    min_len = models.IntegerField(verbose_name='min_len')
    max_len = models.IntegerField(verbose_name='max_len')
    assembly_time = models.IntegerField(verbose_name='assembly_time')
    verification_time = models.IntegerField(verbose_name='verification_time ')
