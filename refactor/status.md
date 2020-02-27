# Problems

## mem leaks 


### 1.3.15_dev_double

```
./bin/geneid -G -3 -P param/dictyostelium.param  ../genomes_4geneidc_t/d.disc/Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa  > dev1.3.15_double__memsan_dict_param_dict_whole_genome_.gff
```

```
=================================================================
==29866==ERROR: LeakSanitizer: detected memory leaks

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5955 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:154:26
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d545b in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:96:27
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5190 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:47:26
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5c4d in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:203:27
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5471 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:103:34
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d51d0 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:61:34
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d549c in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:110:34
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d51a6 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:54:34
    #2 0x55e49d8d79dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:483:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5c63 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:210:34
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5996 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:168:34
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d5c8f in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:217:34
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x55e49d874679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x55e49d8d596b in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:161:34
    #2 0x55e49d8d8081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:532:9
    #3 0x55e49d8aa5e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7feff154a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

SUMMARY: AddressSanitizer: 350496 byte(s) leaked in 396 allocation(s).
```
