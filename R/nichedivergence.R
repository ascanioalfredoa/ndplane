C <-
    A[A$x %in% B$x, ] %>%
    left_join(B[B$x %in% A$x, ], by = "x", suffix = c("_a", "_b")) %>%
    pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") %>%
    group_by(x) %>%
    summarise(y = min(y))

NicheDiss <- 1 - (trapz(C$x, C$y)/trapz(A$x, A$y) + trapz(C$x, C$y)/trapz(B$x, B$y))/2

RNO <- 1 - (min(c(max(A$x), max(B$x))) - max(c(min(A$x), min(B$x))))/(max(c(max(A$x), max(B$x))) - min(c(min(A$x), min(B$x))))
if(RNO > 1) RNO <- 1
