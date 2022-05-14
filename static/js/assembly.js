new Vue({
    el: '#app',
    delimiters: ['{[', ']}'],
    data() {
        return {
            output1: false,
            output2: false,

            showTable: true,
            loading: false,
            verificationLoading: false,
            another_res: false,

            activeIndex: 'second',
            activeName: 'input',

            dynamicValidateForm: {
                email: null,
                geneLen: 0,
                result: 'res1',
                minLen: 20,
                maxLen: 30,
                pools: 1,
                temperature: 54,
                oligoConc: 100,
                primerConc: 100,

                // concentrations: 1,

                resultType: 'Gap',
                verification: 'No',

                tableData: [
                    {
                        name: 'K',
                        data: 50,
                        ion: 'Na<sup>+</sup> / K<sup>+</sup> (1-1000mM)',
                        unit: 'mM',
                        info: '',
                    }, {
                        name: 'Mg',
                        data: 8,
                        ion: 'Mg<sup>2+</sup> (1-1000mM)',
                        unit: 'mM',
                        info: 'The Mg-dNTP association constant is sufficiently large that the free Magnesium ion concentration can be approximated simply by the difference between the total magnesium concentration and the concentration of dNTPs.',
                    }, {
                        name: 'dNTPs',
                        data: 4,
                        ion: 'dNTPs (1-1000mM)',
                        unit: 'mM',
                        info: '',
                    }, {
                        name: 'Tris',
                        data: 10,
                        ion: 'Tris<sup>2+</sup> (1-1000mM)',
                        unit: 'mM',
                        info: '',
                    }, {
                        name: 'oligo',
                        data: 10,
                        ion: 'oligo (0.1-1000nM)',
                        unit: 'nM',
                        info: '',
                    }, {
                        name: 'primer',
                        data: 400,
                        ion: 'primer (0-1000nM)',
                        unit: 'nM',
                        info: '',
                    }],
                geneName: null,
                geneDesc: null,

                gene: null,
            },

            dataColumn: [
                {
                    name: 'index', width: 95,
                },{
                    name: "oligo(5'->3')", width: 633,
                },{
                    name: 'Tm', width: 70,
                },{
                    name: 'overlap', width: 80,
                },{
                    name: 'length', width: 75,
                },
            ],

            arr: [],
            anotherArr: [],
        }
    },

    created() {
        // {#this.calGeneLength();#}
        // this.dynamicValidateForm.email = '758168660@qq.com';
        // this.dynamicValidateForm.geneDesc = 'description';
        // this.dynamicValidateForm.geneName = 'name';
        // this.dynamicValidateForm.gene = 'taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaataataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata';
        // // this.dynamicValidateForm.gene = 'taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggtaagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacataagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccggttgcgcagcctgaataataacgctaagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagctgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaactgggcctttcgttttat';
        // this.dynamicValidateForm.geneLen = this.dynamicValidateForm.gene.length;
    },

    methods: {
        handleSelect(key, keyPath) {
            // var that = this;
            this.activeIndex = key;
            // console.log(this.activeIndex);
            // this.$emit.activeIndex = key;
        },

        test(resArr) {
            // view tail
            source = resArr[resArr.length - 1];
            // console.log(source)
            var tail = source["tail_reverse"];
            if (tail.length > 0) {
                // console.log(tail);
                replaceStr = '<span style="color:red;">' + tail + '</span>';
                index = source["info"].length - 3;
                str = source["info"][index][1];
                str = str.replace(tail, replaceStr);
                // str = str.substring(0, str.lastIndexOf(tail)) + replaceStr

                source["info"][index][1] = str;
            }

            var above_tail = source["above_tail"];
            if (above_tail.length > 0) {
                // console.log(above_tail);
                replaceStr = '<span style="color:red;">' + above_tail + '</span>';
                index = source["info"].length - 4;
                str = source["info"][index][1];
                str = str.replace(above_tail, replaceStr);
                // str = str.substring(0, str.lastIndexOf(tail)) + replaceStr

                source["info"][index][1] = str;
            }

            return resArr;
        },

        dateFormat(fmt, date) {
            let ret;
            const opt = {
                "Y+": date.getFullYear().toString(),        // 年
                "m+": (date.getMonth() + 1).toString(),     // 月
                "d+": date.getDate().toString(),            // 日
                "H+": date.getHours().toString(),           // 时
                "M+": date.getMinutes().toString(),         // 分
                "S+": date.getSeconds().toString()          // 秒
                // 有其他格式化字符需求可以继续添加，必须转化成字符串
            };

            for (let k in opt) {
                ret = new RegExp("(" + k + ")").exec(fmt);
                if (ret) {
                    fmt = fmt.replace(ret[1], (ret[1].length === 1) ? (opt[k]) : (opt[k].padStart(ret[1].length, "0")))
                }
            }

            return fmt;
        },

        download(res) {
            var that = this;
            var input = that.arr;

            if (res == 'res2') {
                input = that.anotherArr;
            }
            // console.log(input);
            // {#console.log(this.arr);#}
            axios.post("/download/", input, {responseType: 'blob'}).then(function (response) {
                // {#console.log(response);#}
                var time = that.dateFormat("YYYY-mm-dd-HH:MM", new Date());
                // {#console.log(time);#}
                const blob = new Blob([response.data]);//处理文档流
                const fileName = that.dynamicValidateForm.geneName + "_" + time + '.xlsx';
                const elink = document.createElement('a');
                elink.download = fileName;
                elink.style.display = 'none';
                elink.href = URL.createObjectURL(blob);
                document.body.appendChild(elink);
                elink.click();
                URL.revokeObjectURL(elink.href); // 释放URL 对象
                document.body.removeChild(elink);

            }).catch(function (error) {
                alert(error);
                // console.log(error);
            });
        },

        testDemo() {
            this.dynamicValidateForm.email = '758168660@qq.com';
            this.dynamicValidateForm.geneDesc = 'description';
            this.dynamicValidateForm.geneName = 'name';
            this.dynamicValidateForm.gene = 'taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaataataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata';
            // this.dynamicValidateForm.gene = 'taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggtaagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacataagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagcttatcgataccggttgcgcagcctgaataataacgctaagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactattttacctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgatatcaagctgcgcaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggcagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaactgggcctttcgttttat';
            this.dynamicValidateForm.geneLen = this.dynamicValidateForm.gene.length;
        },

        lastResult() {
            this.activeName = 'output1';
        },

        refreshTable() {
            let that = this;
            that.showTable = false;
            that.$nextTick(() => that.showTable = true);
        },

        resetForm(formName) {
            this.$refs[formName].resetFields();
        },

        testData(formData) {
            // email
            if (formData.email == null) {
                this.$message({
                        message: 'Please enter your email address',
                        type: "error"
                    }
                );
                return false;
            }

            // gene
            if (formData.gene == null) {
                this.$message({
                    message: 'please enter the gene seq',
                    type: "error"
                });
                return false;
            }

            // overlap
            if (formData.minLen <= 0 || formData.minLen >= formData.maxLen || formData.maxLen >= formData.geneLen) {
                this.$message({
                    message: 'overlap length error ',
                    type: "error"
                });
                return false;
            }

            // Parameters
            for (let i = 0; i < 4; i++) {
                if (formData.tableData[i].data < 1 || formData.tableData[i].data > 1000) {
                    this.$message({
                        message: formData.tableData[i].name + ' ion concentration',
                        type: "error"
                    });
                    return false;
                }
            }
            if (formData.tableData[4].data < 0.1 || formData.tableData[4].data > 1000) {
                this.$message({
                    message: formData.tableData[4].name + ' ion concentration',
                    type: "error"
                });
                return false;
            }
            if (formData.tableData[5].data < 0 || formData.tableData[5].data > 1000) {
                this.$message({
                    message: formData.tableData[5].name + ' ion concentration',
                    type: "error"
                });
                return false;
            }
            if (formData.tableData[1].data - formData.tableData[2].data <= 0) {
                this.$message({
                    message: 'Magnesium ion concentration',
                    type: "error"
                });
                return false;
            }

            // pools
            if (formData.pools <= 0) {
                this.$message({
                    message: 'pools error ',
                    type: "error"
                });
                return false;
            }

            return true
        },

        submitForm(formData) {

            if(!this.testData(formData)) {
                return;
            }

            var that = this;
            // {#console.log(formData);#}
            that.loading = true;
            that.output1 = false;
            that.output2 = false;
            that.arr = [];
            that.anotherArr = [];
            formData.result = 'res1';

            if (that.dynamicValidateForm.pools == 1) {
                axios.post("/assembly/", formData).then(function (response) {
                    that.arr = response.data.arr;
                    // console.log(that.arr)
                    that.arr = that.test(that.arr);

                    that.activeName = 'output1';
                    that.output1 = true;
                    that.loading = false;
                }).catch(function (error) {
                    // console.log(error);
                    alert(error);

                    that.loading = false;
                });

            } else {
                axios.post("/assemblyPools/", formData).then(function (response) {
                    that.arr = response.data.arr;
                    // console.log(that.arr)
                    that.arr = that.test(that.arr);
                    // console.log(response.data.arr);

                    that.activeName = 'output1';
                    that.output1 = true;
                    that.loading = false;

                }).catch(function (error) {
                    // console.log(error);
                    alert(error);
                    that.loading = false;
                });
            }
        },

        validation(index, nextCal, tem, oligoConc, primerConc, res) {
            var that = this;
            var info = [nextCal, tem, oligoConc, primerConc]

            that.verificationLoading = true;
            that.loading = true;
            axios.post("/analysis/", info).then(function (response) {
                if(res == 'res1') {
                    that.arr[index]["analyInfo"] = response.data.analyInfo;
                } else if(res == 'res2') {
                    that.anotherArr[index]["analyInfo"] = response.data.analyInfo;
                }

                that.dynamicValidateForm.pools = 1;
                that.$nextTick(() => that.dynamicValidateForm.pools = 2);

                that.verificationLoading = false;
                that.loading = false;

            }).catch(function (error) {
                alert(error);
                that.verificationLoading = false;
                that.loading = false;
                // console.log(error);
            });

        },

        anotherResult(formData) {
            if (this.anotherArr.length !== 0) {
                this.activeName = 'output2';

            } else {
                var that = this;
                formData.result = 'res2';
                that.another_res = true;

                if (that.dynamicValidateForm.pools == 1) {
                    axios.post("/assembly/", formData).then(function (response) {
                        that.anotherArr = response.data.arr;
                        that.anotherArr = that.test(that.anotherArr)

                        // console.log(that.anotherArr)
                        that.activeName = 'output2';
                        that.output2 = true;
                        that.loading = false;
                        that.another_res = false;
                    }).catch(function (error) {
                        // console.log(error);
                        alert(error);
                        that.loading = false;
                        that.another_res = false;
                    });
                } else {
                    axios.post("/assemblyPools/", formData).then(function (response) {
                        that.anotherArr = response.data.arr;
                        // console.log(that.anotherArr);
                        that.anotherArr = that.test(that.anotherArr);
                        that.activeName = 'output2';
                        that.output2 = true;
                        that.loading = false;
                        that.another_res = false;

                    }).catch(function (error) {
                        // console.log(error);
                        alert(error);
                        that.loading = false;
                    });
                }
            }
        },

        calGeneLength() {
            this.dynamicValidateForm.geneLen = this.dynamicValidateForm.gene.length;
        },
    },


});
