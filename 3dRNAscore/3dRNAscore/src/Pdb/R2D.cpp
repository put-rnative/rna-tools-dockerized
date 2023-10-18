#include "R2D.h"

using namespace ss;

RNA_Module::RNA_Module(int type, int len, string seq, string ss, int *num, int flag) {
	type_ = type;
	len_ = len;
	seq_ = seq;
	ss_ = ss;
	num_ = new int[len];
	for (int i = 0; i < len; i++) {
		num_[i] = num[i];
	}
	flag_ = flag;
	son_ = NULL;
	brother_ = NULL;
}

RNA_Module::~RNA_Module() {
	delete [] num_;
}

R2D::R2D(string seq, string second_struct, int view) {
	// check sequence
	int len_seq = 0;
	for (int i = 0; i < seq.size(); i++) {
		if (seq[i] == 'A' || seq[i] == 'a' || 
		seq[i] == 'U' || seq[i] == 'u' || 
		seq[i] == 'G' || seq[i] == 'g' || 
		seq[i] == 'C' || seq[i] == 'c') {
			len_seq++;
		} else {
			cerr << "RNA sequence has no residue named '" << seq[i] << "'" << endl;
			exit(1);
		}
	}
	seq_ = seq;

	// check secondary structure
	int len_second_struct = 0;
	ss_ = "";
	string second_struct1 = ""; // secondary structure which has no '[' and ']';
	string second_struct2 = ""; // secondary structure which has no '(' and ')';
	for (int i = 0; i < second_struct.size(); i++) {
		if (second_struct[i] == '.') {
			len_second_struct++;
			ss_ += '.';
			second_struct1 += '.';
			second_struct2 += '.';
		} else if (second_struct[i] == '(' || second_struct[i] == ')') {
			len_second_struct++;
			ss_ += second_struct[i];
			second_struct1 += second_struct[i];
			second_struct2 += '.';
		} else if (second_struct[i] == '[' || second_struct[i] == ']') {
			len_second_struct++;
			ss_ += '.';
			second_struct1 += '.';
			if (second_struct[i] == '[') { 
				second_struct2 += '('; // transform '[' to '(' 
			} else {
				second_struct2 += ')'; // transform ']' to ')'
			}
		} else if (second_struct[i] == '&') {
		} else {
			cerr << "Please input a legal dot-bracket form of secondary structure!" << endl;
			exit(1);
		}
	}
	if (len_seq != len_second_struct) {
		cerr << "The length of the sequence and the length of the secondary structure must be equal!" << endl;
		exit(1);
	}
	len_ = len_seq;

	view_ = view;

	// set 2d structure tree
	head_ = setTree(seq_, second_struct1);
	pseudo_head_ = setTree(seq_, second_struct2);
}

RNA_Module* R2D::setTree(string seq, string second_struct) {
	int len = seq.size();
	vector<Nucleotide> left_stack;
	vector<RNA_Module *> module_stack;
	RNA_Module *head;
	for (int i = 0; i < second_struct.size(); i++) {
		Nucleotide nuc(seq[i], second_struct[i], i);
		left_stack.push_back(nuc);
		if (i == 0) {
			if (nuc.ss_ == ')') {
				cerr << "Wrong secondary structure!" << endl;
				exit(1);
			}
		} else {
			if (nuc.ss_ == ')') {
				if ((i + 1 < len && second_struct[i + 1] != ')') || i + 1 == len) {
					vector<Nucleotide> right_stack;
					int left_bracket_number = 0;
					int right_bracket_number = 0;
					while (1) {
						Nucleotide temp_nuc = left_stack.back();
						right_stack.push_back(temp_nuc);
						left_stack.pop_back();
						if (temp_nuc.ss_ == '(') {
							left_bracket_number++;
						} else if (temp_nuc.ss_ == ')') {
							right_bracket_number++;
						}
						if (temp_nuc.ss_ == '(' && ((left_stack.back().ss_ != '(') || left_stack.empty())) {
							vector<Nucleotide> temp_stack;
							int temp_left_number = 0;
							int temp_right_number = 0;
							int p1 = 0, p2, p3, p4;
							int loop_type = 0;
							int	loop_length = 0;
							int j = 0;
							while (1) {
								Nucleotide temp_nuc2 = right_stack.back();
								temp_stack.push_back(temp_nuc2);
								right_stack.pop_back();
								if (temp_nuc2.ss_ == '(') {
									temp_left_number++;
									if (temp_left_number == left_bracket_number) {
										p2 = j;
									}
								} else if (temp_nuc2.ss_ == ')') {
									temp_right_number++;
									if (temp_right_number == 1) {
										p3 = j;
									} else if (temp_right_number == left_bracket_number) {
										p4 = j;
									}
								} else if (temp_nuc2.ss_ == '[') {
									loop_type++;
								} else if (temp_nuc2.ss_ == ']') {
								} else if (temp_nuc2.ss_ == '.') {
									loop_length++;
								}
								if (temp_left_number == temp_right_number) {
									int temp_length;
									string temp_seq;
									string temp_ss;
									int *temp_num;

									if (loop_length != 0) {
										temp_length = p3 - p2 + 1;
										temp_num = new int[temp_length];
										for (int k = 0; k < temp_length; k++) {
											temp_seq += temp_stack[k + p2].seq_;
											if (temp_stack[k + p2].ss_ == '[') {
												temp_ss += '(';
											} else if (temp_stack[k + p2].ss_ == ']') {
												temp_ss += ')';
											} else {
												temp_ss += temp_stack[k + p2].ss_;
											}
											temp_num[k] = temp_stack[k + p2].num_;
										}
										RNA_Module *temp_loop = new RNA_Module(loop_type + 1, temp_length, temp_seq, temp_ss, temp_num);
										if (loop_type > 0) {
											RNA_Module *temp_module = module_stack.back();
											temp_loop->son_ = temp_module;
											module_stack.pop_back();
											for (int k = 0; k < loop_type - 1; k++) {
												RNA_Module *new_temp_module = module_stack.back();
												temp_module->brother_ = new_temp_module;
												module_stack.pop_back();
												temp_module = new_temp_module;
											}
											module_stack.push_back(temp_loop);
										} else {
											module_stack.push_back(temp_loop);
										}
										head = temp_loop;
									}

									temp_length = temp_left_number + temp_right_number;
									temp_num = new int[temp_length];
									temp_seq = "";
									temp_ss = "";
									for (int k = 0; k < temp_left_number; k++) {
										temp_seq += temp_stack[k].seq_;
										temp_ss += '(';
										temp_num[k] = temp_stack[k].num_;
									}
									for (int k = temp_left_number; k < temp_length; k++) {
										temp_seq += temp_stack[k + p3 - p2 - 1].seq_;
										temp_ss += ')';
										temp_num[k] = temp_stack[k + p3 - p2 - 1].num_;
									}
									RNA_Module *temp_helix = new RNA_Module(0, temp_length, temp_seq, temp_ss, temp_num);
									if (!module_stack.empty()) {
										RNA_Module *temp_module = module_stack.back();
										temp_helix->son_ = temp_module;
										module_stack.pop_back();
									}
									module_stack.push_back(temp_helix);
									head = temp_helix;

									Nucleotide temp_nuc1(temp_nuc2.seq_, ']', temp_nuc2.num_);
									right_stack.push_back(temp_nuc1);
									Nucleotide temp_nuc2(temp_nuc.seq_, '[', temp_nuc.num_);
									right_stack.push_back(temp_nuc2);
									break;
								}
								j++;
							}
						}
						if (left_bracket_number == right_bracket_number) {
							Nucleotide temp_nuc1(temp_nuc.seq_, '[', temp_nuc.num_);
							left_stack.push_back(temp_nuc1);
							Nucleotide temp_nuc2(nuc.seq_, ']', nuc.num_);
							left_stack.push_back(temp_nuc2);
							break;
						}
						if (left_stack.empty() && left_bracket_number != right_bracket_number) {
							cerr << "Wrong secondary structure!" << endl;
							exit(1);
						}
					}
				}
			}
		}
	}
	if (!left_stack.empty()) {
		int loop_type = 0;
		int temp_length = left_stack.size();
		string temp_seq;
		string temp_ss;
		int *temp_num = new int[temp_length];
		for (int i = 0; i < temp_length; i++) {
			if (left_stack[i].ss_ == '(' || left_stack[i].ss_ == ')') {
				cerr << "Wrong secondary structure!" << endl;
				exit(1);
			} else if (left_stack[i].ss_ == '[') {
				loop_type++;
				temp_ss += '(';
			} else if (left_stack[i].ss_ == ']') {
				temp_ss += ')';
			} else if (left_stack[i].ss_ == '.') {
				temp_ss += '.';
			}
			temp_seq += left_stack[i].seq_;
			temp_num[i] = left_stack[i].num_;
		}
		RNA_Module *temp_loop = new RNA_Module(-loop_type - 1, temp_length, temp_seq, temp_ss, temp_num);
		if (temp_ss == "()") {
			return head;
		}
		if (!module_stack.empty()) {
			RNA_Module *temp_module = module_stack.back();
			temp_loop->son_ = temp_module;
			module_stack.pop_back();
			for (int i = 0; i < loop_type - 1; i++) {
				RNA_Module *new_temp_module = module_stack.back();
				temp_module->brother_ = new_temp_module;
				module_stack.pop_back();
				temp_module = new_temp_module;
			}
		}
		module_stack.push_back(temp_loop);
		head = temp_loop;
	}
	return head;
}

R2D::~R2D() {
	delTree(head_);
	delTree(pseudo_head_);
}

void R2D::delTree(RNA_Module *head) {
	if (head == NULL) {
		return;
	}
	delTree(head->son_);
	delTree(head->brother_);
	delete head;
}

void R2D::print() {
	cout << "================= Secondary structure tree =====================" << endl;
	printTree(head_);
	cout << endl;
	cout << "================= Pseudoknot tree ======================" << endl;
	printTree(pseudo_head_);
	cout << endl;
}

void R2D::printTree(RNA_Module *head) {
	if (head == NULL) {
		return;
	}

	cout << head->seq_ << endl;
	cout << head->ss_ << endl;
	for (int i = 0; i < head->len_; i++) {
		cout << head->num_[i] << ' ';
	}
	cout << endl;
	cout << endl;

	printTree(head->son_);
	printTree(head->brother_);
}

