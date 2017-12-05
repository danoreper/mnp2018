mysql -u root -p -e '

select fm.*, strain_name 

from founder_meta fm 

inner join strain on fm.strain_id=strain.strain_id 
limit 10000000000;' paul_1410 > out2.txt

