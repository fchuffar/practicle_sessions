##### 1. Creating account on perseus
##### https://perseus.ujf-grenoble.fr/

##### 2. Fine tunning of ssh
# Generate keys, with passphrase for cluster access
ssh-keygen
# Generate keys, without passphrase for accessing current jobs
ssh-keygen -f ~/.ssh/id_rsa_open
# Set ssh configuration
# use access-ciment machine to access to ciment clusters
echo "
Host *
   StrictHostKeyChecking no
   HashKnownHosts no
 Host rotule
   User fchdemo
   IdentityFile  ~/.ssh/id_rsa
   Hostname rotule.univ-grenoble-alpes.fr
 Host trinity
   User fchdemo
   IdentityFile  ~/.ssh/id_rsa
   Hostname trinity.univ-grenoble-alpes.fr
Host luke dahu cargo
   User fchdemo
   ProxyCommand ssh fchdemo@access-gricad.univ-grenoble-alpes.fr -W %h:%p
" > .ssh/config 
# Copy keys and config on access-ciment (trinity and rotule)
scp -r .ssh trinity:. 
scp -r .ssh rotule:.

# Connection on access-ciment to authorize keys and copy config to dahu
ssh trinity
cat .ssh/id_rsa.pub > .ssh/authorized_keys
ssh rotule
cat .ssh/id_rsa.pub > .ssh/authorized_keys
scp -r .ssh dahu:.
scp -r .ssh luke:.
scp -r .ssh cargo:.

##### 3. Fine tunning of ~/.profile
# On dahu (ssh dahu)
# Add epimed environment to your PATH
echo 'export PATH="/summer/epistorage/opt/bin:$PATH"' >> ~/.profile
echo 'export PATH="/summer/epistorage/miniconda3/bin:$PATH"' >> ~/.profile
# export OAR_JOB_KEY_FILE (useful without passphrase to automate supervisation)
echo 'export OAR_JOB_KEY_FILE=~/.ssh/id_rsa' >> ~/.profile

##### 4. launch resources
# Interactive connection to ressources
oarsub -I --project epimed
# Launch batch on computing element
oarsub --project epimed "sleep 3600"
# Lauch devel job (30min on a node, i.e., 32 cores)
oarsub --project epimed  -l nodes=1,walltime=00:30:00 -t devel -I
