##### 1. Creating account on perseus
##### https://perseus.ujf-grenoble.fr/

##### 2. Fine tunning of ssh
# Generate keys with passphrase
ssh-keygen
# Set ssh configuration
# use access-ciment machine to access to ciment clusters
echo "
Host *
   StrictHostKeyChecking no
   HashKnownHosts no
Host trinity
   User fchdemo
   hostname trinity.ujf-grenoble.fr
Host rotule
   User fchdemo
   hostname rotule.imag.fr
Host luke dahu
   User fchdemo
   ProxyCommand ssh fchdemo@access-ciment.imag.fr -W %h:%p
" > ~/.ssh/config 
# Copy keys and config on access-ciment (trinity and rotule)
scp -r .ssh fchdemo@access-ciment.imag.fr:.
scp -r .ssh trinity:. 
scp -r .ssh rotule:.

# Connection on access-ciment to authorize keys and copy config to luke
ssh fchdemo@access-ciment.imag.fr
cat .ssh/id_rsa.pub > .ssh/authorized_keys
scp -r .ssh luke:.
scp -r .ssh dahu:.

##### 3. Fine tunning of ~/.profile
# On dahu (ssh dahu)
ssh dahu
# Add epimed environment to your PATH
echo 'export PATH="/summer/epistorage/opt/bin:$PATH"' >> ~/.profile
echo 'export PATH="/summer/epistorage/miniconda3/bin:$PATH"' >> ~/.profile
# export OAR_JOB_KEY_FILE (useful without passphrase to automate supervisation)
echo 'export OAR_JOB_KEY_FILE=~/.ssh/id_rsa' >> ~/.profile

##### 4. launch resources
# overview of available ressources
chandler
# Interactive connection to ressources
oarsub -I --project epimed
# Launch batch on computing element
# oarsub --project epimed "sleep 30"

