from django.urls import path
from app import views

#Template tagging
app_name = 'app'

urlpatterns = [
    path('', views.predict, name='predict')
]