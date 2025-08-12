from django.urls import path
from . import views
from django.conf.urls.static import static
from django.conf import settings
urlpatterns = [
    path('', views.upload_protein, name='upload_protein'),
    # Add other URL patterns as needed
]+static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
