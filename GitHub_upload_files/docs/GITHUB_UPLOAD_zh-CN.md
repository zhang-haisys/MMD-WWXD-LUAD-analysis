# 上传GitHub

1. 解压交付的ZIP；进入MMD_WXD_GitHub_code目录。
2. 在GitHub新建仓库，填写研究对应名称。此代码尚未完成复现验证，
   仓库描述建议使用“Analysis source code and reproducibility documentation”。
3. 保留README、docs、metadata及原始脚本，上传本目录内容。
   使用Git客户端上传可保留.gitignore等隐藏文件，ZIP本身不是代码仓库布局。
4. 不上传config.local.R、workspace内本地数据或results中的敏感运行信息。
5. 作者确认第三方脚本的再分发权限后补齐R/vendor，选择批准的许可证。
6. 补齐并验证后创建明确版本的release；论文填写真实仓库网址/版本或归档DOI。
   不应预先编造网址、DOI或承诺完全可复现。

可选Git命令（在仓库根目录运行；本次未执行）：

```text
git init
git add .
git commit -m "Organize supplied study scripts and document dependencies"
git branch -M main
git remote add origin <真实仓库地址>
git push -u origin main
```

若上传的是本次未完成复现的版本，代码可用性声明应明确其范围：
“The supplied analysis scripts and documentation are available at [verified URL].
Dependencies and limitations of historical reproducibility are documented in the repository.”
不要在尚未补齐辅助脚本时声称“all code required to reproduce the results”。
