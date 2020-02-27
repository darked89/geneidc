# Problems

## mem leaks 


### 1.3.15_dev

```
./bin/geneid -G -3 -P param/dictyostelium.param  ../genomes_4geneidc_t/d.disc/Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa  > dev1.3.15__memsan_dict_param_dict_whole_genome.gff
```

```
=================================================================
==356==ERROR: LeakSanitizer: detected memory leaks

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd955 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:155:26
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd45b in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:97:27
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd190 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:48:26
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Direct leak of 83304 byte(s) in 39 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcdc4d in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:204:27
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd471 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:104:34
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd1d0 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:62:34
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd49c in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:111:34
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 2736 byte(s) in 38 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd1a6 in InsertBeginExon /home/darked/proj/geneidc/./src/SortExons.c:55:34
    #2 0x561ae5dcf9dc in SortExons /home/darked/proj/geneidc/./src/SortExons.c:484:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcdc63 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:211:34
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd996 in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:169:34
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcdc8f in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:218:34
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

Indirect leak of 1584 byte(s) in 22 object(s) allocated from:
    #0 0x561ae5d6c679 in malloc (/home/darked/proj/geneidc/bin/geneid+0xd7679)
    #1 0x561ae5dcd96b in InsertEndExon /home/darked/proj/geneidc/./src/SortExons.c:162:34
    #2 0x561ae5dd0081 in SortExons /home/darked/proj/geneidc/./src/SortExons.c:533:9
    #3 0x561ae5da25e0 in main /home/darked/proj/geneidc/./src/geneid.c:515:17
    #4 0x7f7df880a152 in __libc_start_main (/usr/lib/libc.so.6+0x27152)

SUMMARY: AddressSanitizer: 350496 byte(s) leaked in 396 allocation(s).

```
