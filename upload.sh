echo "uploading files"
rsync -azv -P \
  --exclude=.git \
  --exclude=.idea \
  --exclude=node_modules \
  --exclude=generated-images \
  --whole-file \
  ~/code/repos/smilesDrawer "su0742@horeka.scc.kit.edu:~/code/repos"
echo "uploading files done"