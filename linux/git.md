# git

## 全局忽略文件

```bash
echo "*.aux" >> ~/.gitignore_global
git config --global core.excludesfile ~/.gitignore_global
```