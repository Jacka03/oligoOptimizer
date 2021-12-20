from django.urls import path

from analysis import views

urlpatterns = [
    path('', views.AssemblyView.as_view()),
    path('assembly/', views.AssemblyView.as_view()),
]