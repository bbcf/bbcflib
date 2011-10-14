make clean html
/usr/bin/ssh -t ginger.epfl.ch "sudo /usr/bin/scp -r '"$USER"@"$HOST":"$PWD"/_build/html/*' /srv/bbcflib"
