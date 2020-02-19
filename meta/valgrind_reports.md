# Fixing memory leaks

```
valgrind -v --track-origins=yes --leak-check=yes  ./bin/geneid   -GP param/dictyostelium.param ./samples/dict_1chr.fa

```


geneid v1.4.5 (upstream repo)
LEAK SUMMARY:
==30662==    definitely lost: 64,528 bytes in 14 blocks
==30662==    indirectly lost: 0 bytes in 0 blocks
==30662==      possibly lost: 0 bytes in 0 blocks
==30662==    still reachable: 2,920,082,070 bytes in 806 blocks
==30662==         suppressed: 0 bytes in 0 blocks
==30662== Reachable blocks (those to which a pointer was found) are not shown.
==30662== To see them, rerun with: --leak-check=full --show-leak-kinds=all


geneidc 20200212
LEAK SUMMARY:
==29127==    definitely lost: 8,512 bytes in 4 blocks
==29127==    indirectly lost: 224 bytes in 4 blocks
==29127==      possibly lost: 0 bytes in 0 blocks
==29127==    still reachable: 6,613,923,270 bytes in 806 blocks
==29127==         suppressed: 0 bytes in 0 blocks
==29127== Reachable blocks (those to which a pointer was found) are not shown.
==29127== To see them, rerun with: --leak-check=full --show-leak-kinds=all

