# VS Code

## 某些问题的解决方法 
<!-- Approach for some issues -->

1. 鼠须管输入法输入中文

   ref: https://github.com/rime/squirrel/issues/179#issuecomment-364048969

   添加自定义squirrel.custom.yaml

   ```
   # Squirrel settings
   # encoding: utf-8

   patch:
     app_options:
       com.microsoft.VSCode:
         ascii_mode: false
   ```

2. python + pylance 某些符号不能识别

   ref: https://github.com/microsoft/pylance-release/issues/78#issuecomment-774516948

   ```
   "python.analysis.extraPaths": [
       "vendor/django-allauth",
       "vendor/dj-rest-auth",
   ],
   ```