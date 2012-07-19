#!/usr/bin/perl
while (<>) {
    s/^@/>/;
    print;
    $_=<>;
    print;
    <>;
    <>;
}
