make clean html
ssh -t ginger.epfl.ch "cd /srv/bbcflib; sudo scp -r '"$USER"@"$HOST":"$PWD"/_build/html/*' ."
