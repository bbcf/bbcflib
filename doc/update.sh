make clean html
ssh -t ginger.epfl.ch "sudo /usr/bin/scp -r '"$USER"@"$HOST":"$PWD"/_build/html/*' /srv/bbcflib"
