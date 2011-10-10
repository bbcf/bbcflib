make clean html
sudo chown -R apache:apache _build/html
rsync -av --owner --group _build/html/ root@ginger.epfl.ch:/srv/bbcflib
sudo chown -R YOURUSERNAME:YOURGROUP _build/html
