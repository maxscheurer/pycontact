create table residues(name TEXT, scpolarity TEXT, hbondtype TEXT);
/* unpolar */ 
insert into residues(name, scpolarity, hbondtype) values("gly","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("ala","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("leu","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("ile","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("val","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("phe","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("met","nonpolar","none");
insert into residues(name, scpolarity, hbondtype) values("pro","nonpolar","none");
/* charged */
insert into residues(name, scpolarity, hbondtype) values("asp","negative","acc");
insert into residues(name, scpolarity, hbondtype) values("glu","negative","acc");
insert into residues(name, scpolarity, hbondtype) values("lys","positive","don");
insert into residues(name, scpolarity, hbondtype) values("arg","positive","don");
insert into residues(name, scpolarity, hbondtype) values("hsp","positive","don");
/* polar */
insert into residues(name, scpolarity, hbondtype) values("asn","polar","both");
insert into residues(name, scpolarity, hbondtype) values("gln","polar","both");
/* histidines uncharged */
insert into residues(name, scpolarity, hbondtype) values("his","polar","both");
insert into residues(name, scpolarity, hbondtype) values("hsd","polar","both");
insert into residues(name, scpolarity, hbondtype) values("hse","polar","both");

insert into residues(name, scpolarity, hbondtype) values("cys","polar","donor");
insert into residues(name, scpolarity, hbondtype) values("ser","polar","both");
insert into residues(name, scpolarity, hbondtype) values("thr","polar","both");
/* special */
insert into residues(name, scpolarity, hbondtype) values("tyr","nonpolar","donor");
insert into residues(name, scpolarity, hbondtype) values("trp","nonpolar","donor");
