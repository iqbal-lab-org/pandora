echo "##############################################################"
echo "START: $(date)"
branch=$(cat /nfs/research1/zi/rmcolq/git/pandora/.git/HEAD | awk '{print $2}')
echo "git branch: $branch"
commit=$(cat /nfs/research1/zi/rmcolq/git/pandora/.git/$branch)
echo "git commit: $commit"
echo "python running: $(which python)"
echo "##############################################################"
