var GridBuilder = (function () {
    "use strict";
    var mIsFirstCall = true,
        mSelf = null,
        mCurrentPath = [],
        mPathTable = [],
        mCellMap = {},
        mTopSequence = "",
        mSideSequence = "",
        mDomGridTable = null,
        mDomAlignmentTable = null,
        mDomContainer = null,
        mDomResultContainer = null,
        mGapSymbol = "-",
        mIsCustomPathMode = false,
        mMatchScore = 0,
        mMismatchScore = 0,
        mGapScore = 0;

    function onCellClicked(dom, x, y) {
        x = parseInt(x, 10);
        y = parseInt(y, 10);

        var lastElement = null;
        if (mCurrentPath !== null && mCurrentPath.length !== 0) {
            lastElement = mCurrentPath[mCurrentPath.length - 1];

            if (dom.hasClass('in-path')) {
                if (mCurrentPath.length === 1) {
                    mCurrentPath[0].dom.removeClass('in-path');
                    mCurrentPath[0].dom.removeClass('is-last');
                    mCurrentPath[0].dom.removeAttr('data-index');
                    mCurrentPath = [];
                    onPathUpdate();
                    return true;
                }

                var indexInPath = parseInt(dom.attr('data-index'), 10);
                for (var i = indexInPath + 1; i < mCurrentPath.length; i++) {
                    mCurrentPath[i].dom.removeClass('in-path');
                    mCurrentPath[i].dom.removeClass('is-last');
                    mCurrentPath[i].dom.removeAttr('data-index');
                }
                mCurrentPath.splice(indexInPath + 1, mCurrentPath.length - indexInPath);
                mCurrentPath[mCurrentPath.length - 1].dom.addClass('is-last');
                onPathUpdate();
                return true;
            }

            if (lastElement.x < x || lastElement.y < y) {
                return false;
            }

            if (x - lastElement.x < -1 || y - lastElement.y < -1) {
                return false;
            }
        }

        dom.attr('data-index', mCurrentPath.length);

        mCurrentPath.push({
            'idx': mCurrentPath.length,
            'x': x,
            'y': y,
            'dom': dom,
            'previous': lastElement
        });

        if (lastElement) {
            lastElement.dom.removeClass('is-last');
        }

        dom.addClass('is-last');
        dom.addClass('in-path');
        onPathUpdate();
        return true;
    }

    function onPathUpdate() {
        var alignedTopSeq = '';
        var alignedSideSeq = '';

        $('th').removeClass('included');

        for (var i = mCurrentPath.length - 1; i >= 0; i--) {
            var currentCell = mCurrentPath[i];
            var nextCell = (i > 0) ? mCurrentPath[i - 1] : null;
            var topChar = mTopSequence[currentCell.x - 1];
            var sideChar = mSideSequence[currentCell.y - 1];

            if (!nextCell) {
                continue;
            }

            if (topChar) {
                if (currentCell.x != nextCell.x) {
                    $('#top_seq_' + (currentCell.x - 1)).addClass('included');
                }
            }

            if (sideChar) {
                if (currentCell.y != nextCell.y) {
                    $('#side_seq_' + (currentCell.y - 1)).addClass('included');
                }
            }

            if (nextCell.x - currentCell.x > 0 && nextCell.y - currentCell.y > 0) {
                alignedTopSeq += topChar;
                alignedSideSeq += sideChar;
                continue;
            }

            if (nextCell.x - currentCell.x > 0) {
                sideChar = mGapSymbol;
            }

            if (nextCell.y - currentCell.y > 0) {
                topChar = mGapSymbol;
            }

            alignedTopSeq += topChar;
            alignedSideSeq += sideChar;
        }

        $('#alignment').remove();

        var $table = $('<table />').attr('id', 'alignment');
        mDomAlignmentTable = $table;

        var score = 0;
        var $tr = $('<tr />');
        for (var idxTop in alignedTopSeq) {
            var c1 = alignedTopSeq[idxTop];
            var c2 = alignedSideSeq[idxTop];

            if (c1 === mGapSymbol || c2 === mGapSymbol) {
                score += mGapScore;
            } else if (c1 === c2) {
                score += mMatchScore;
            } else {
                score += mMismatchScore;
            }
            $tr.append($('<td />').html(c1));
        }
        $table.append($tr);

        $tr = $('<tr />');
        for (var idxSide in alignedSideSeq) {
            $tr.append($('<td />').html(alignedSideSeq[idxSide]));
        }
        $table.append($tr);

        $tr = $('<tr />');
        $tr.append($('<td colspan="1500" class="score" />').html("Score = " + score));
        $table.append($tr);

        // mDomResultContainer.append($table);

        // Update result with aligned sequences
        $('#alignedSeq1').text(alignedTopSeq);
        $('#alignedSeq2').text(alignedSideSeq);
        $('#maxScore').text("Maximum Score: " + score);
    }

    function displayTooltip(text, x, y) {
        if ($('#tooltip').length === 0) {
            $('body').prepend($('<div />').attr('id', 'tooltip'));
        }
        var tt = $('#tooltip').html("");
        var tooltipHeight = 30;

        var xBorder = x + tt.width() + 30;
        if (xBorder > $(window).width()) x -= (xBorder - $(window).width());

        var yBorder = y + tt.height() + 30;
        if (yBorder > $(window).height()) y -= (tooltipHeight * 2);

        tt.append(text);
        tt.css('left', x);
        tt.css('top', y);
        tt.css('display', 'block');
    }

    function hideTooltip() {
        $('#tooltip').css('display', 'none');
    }

    function showTooltip(x, y) {
        var targetCell = mCellMap[x + "_" + y];
        var $table = $("<table />");

        var $tr = $("<tr />");
        $tr.append(
            $("<td />").html("<b><u>Score from Diagonal cell</u></b> <br> " + targetCell.diagonalScoreText)
        ).append(
            $("<td />").html("<b><u>Score from Upper cell</u></b> <br> " + targetCell.upScoreText)
        );
        $table.append($tr);

        $tr = $("<tr />");
        $tr.append(
            $("<td />").html("<b><u>Score from Side cell</u></b> <br> " + targetCell.sideScoreText)
        ).append(
            $("<td />").html("Winning (max) score is " + targetCell.winningScore)
        );
        $table.append($tr);

        $('#' + (x - 1) + '_' + (y - 1)).addClass('highlight');
        $('#' + (x - 0) + '_' + (y - 1)).addClass('highlight');
        $('#' + (x - 1) + '_' + (y - 0)).addClass('highlight');

        var targetDom = $('#' + x + '_' + y);
        var pos = targetDom.offset();
        targetDom.addClass('highlight-main');
        displayTooltip($table, pos.left + targetDom.width() + 10, pos.top - targetDom.height() / 2);
    }

    function getCssClassesFromDirection(directions) {
        var cssClasses = "";

        if (!Array.isArray(directions)) {
            return cssClasses;
        }

        cssClasses = directions.join(' ');

        return cssClasses;
    }

    function constructNRow(n) {
        var $table = $('#grid');
        var charIndex = parseInt(n, 10) - 1;
        var $tr = $('<tr />');
        var $th = null;

        if (charIndex >= 0) {
            $th = $('<th />')
                .addClass("seq-header")
                .addClass("side-header")
                .attr('id', 'side_seq_' + charIndex)
                .html(mSideSequence[charIndex]);
            $tr.append($th);
        } else {
            $th = $('<th />');
            $tr.append($th);
        }

        var $td = $('<td />')
            .html(mCellMap[0 + "_" + n].winningScore)
            .attr('data-x', 0)
            .attr('data-y', n)
            .attr('id', 0 + "_" + n);
        $tr.append($td);

        for (var idx in mTopSequence) {
            idx = parseInt(idx, 10);
            var dataPointIndex = (idx + 1) + '_' + (charIndex + 1);

            var cssClasses = "";
            if (n > 0) {
                cssClasses = getCssClassesFromDirection(mCellMap[(idx + 1) + "_" + (charIndex + 1)].direction);
            }

            $td = $('<td />')
                .addClass(cssClasses)
                .html(mCellMap[dataPointIndex].winningScore)
                .attr('data-x', (idx + 1))
                .attr('data-y', (charIndex + 1))
                .attr('data-dg', mCellMap[dataPointIndex].diagonalScoreText)
                .attr('data-up', mCellMap[dataPointIndex].upScoreText)
                .attr('data-sd', mCellMap[dataPointIndex].sideScoreText)
                .attr('id', dataPointIndex);
            $tr.append($td);
        }

        $table.append($tr);
        mDomContainer.append($table);
    }

    function constructGrid() {
        $('#alignment').remove();
        $('#grid').remove();
        var $table = $('<table />').attr('id', 'grid');
        mDomGridTable = $table;
        mDomContainer.append($table);

        var $tr = $('<tr />');

        var $th = $('<th />');
        $tr.append($th);

        $th = $('<th />');
        $tr.append($th);

        for (var idx in mTopSequence) {
            $th = $('<th />');
            $th.attr('id', 'top_seq_' + idx);
            $th.addClass("seq-header");
            $th.addClass("top-header");
            $th.html(mTopSequence[idx]);
            $tr.append($th);
        }

        $table.append($tr);

        for (var i = 0; i < mSideSequence.length + 1; i++) {
            constructNRow(i);
        }

        $('#grid td').click(function() {
            var self = $(this);
            onCellClicked(
                self,
                self.attr('data-x'),
                self.attr('data-y')
            );
        });

        $('#grid td').hover(function() {
            if (mIsCustomPathMode) {
                return;
            }

            var self = $(this);
            var x = self.attr('data-x');
            var y = self.attr('data-y');

            if (x < 1 || y < 1) {
                return;
            }
            $("#side_seq_" + (y - 1)).addClass('highlight');
            $("#top_seq_" + (x - 1)).addClass('highlight');

            showTooltip(x, y);

        }, function() {
            $(".seq-header").removeClass('highlight');
            $('#grid td').removeClass('highlight');
            $('#grid td').removeClass('highlight-main');
            hideTooltip();
        });

        $('#grid th').hover(function() {
            var self = $(this);
            if (!self.hasClass("seq-header")) {
                return;
            }

            var pos = self.offset();
            var topMargin = self.hasClass("side-header") ? self.height() / 4 : self.height() + 4;
            var leftMargin = self.hasClass("side-header") ? self.width() + 4 : 0;
            var text = self.hasClass("included") ? "Included In Alignment" : "Not Included In Alignment";

            displayTooltip(text, pos.left + leftMargin, pos.top + topMargin);

        }, function() {
            hideTooltip();
        });
    }

    mSelf = {
        highlightOptimal: function() {
            mIsCustomPathMode = false;
            var width = mTopSequence.length + 1;
            var height = mSideSequence.length + 1;

            var maxX = 0, maxY = 0, maxScore = 0;

            for (var i = 0; i < width; i++) {
                for (var j = 0; j < height; j++) {
                    if (mPathTable[i][j] > maxScore) {
                        maxScore = mPathTable[i][j];
                        maxX = i;
                        maxY = j;
                    }
                }
            }

            var alignedTopSeq = '';
            var alignedSideSeq = '';

            var currentX = maxX;
            var currentY = maxY;
            while (currentX > 0 && currentY > 0 && mPathTable[currentX][currentY] > 0) {
                var currentCell = mCellMap[currentX + '_' + currentY];
                var currentDom = $('#' + currentX + '_' + currentY);

                alignedTopSeq = mTopSequence[currentX - 1] + alignedTopSeq;
                alignedSideSeq = mSideSequence[currentY - 1] + alignedSideSeq;

                currentDom.click();

                var direction = null;
                if (currentCell.direction) {
                    direction = currentCell.direction[currentCell.direction.length - 1];
                }

                switch (direction) {
                    case 's':
                        currentX--;
                        alignedSideSeq = mGapSymbol + alignedSideSeq;
                        break;
                    case 'u':
                        currentY--;
                        alignedTopSeq = mGapSymbol + alignedTopSeq;
                        break;
                    case 'd':
                        currentX--;
                        currentY--;
                        break;
                    default:
                        break;
                }
            }

            $('#alignedSeq1').text(alignedTopSeq);
            $('#alignedSeq2').text(alignedSideSeq);
            $('#maxScore').text("Maximum Score: " + maxScore);
        },

        startCustomPath: function() {
            this.rebuildTable(mDomContainer, mDomResultContainer, mMatchScore, mMismatchScore, mGapScore, mSideSequence, mTopSequence);
            mIsCustomPathMode = true;
        },

        rebuildTable: function(domContainer, resultContainer, matchScore, mismatchScore, gapScore, seqSide, seqTop) {
            if (mIsFirstCall) {
                $(window).mousemove(function(e) {
                    window.mouseXPos = e.pageX;
                    window.mouseYPos = e.pageY;
                });
                mIsFirstCall = false;
            }

            seqTop = seqTop.toUpperCase();
            seqSide = seqSide.toUpperCase();
            mCurrentPath = [];
            mDomContainer = domContainer;
            mDomResultContainer = resultContainer;
            mTopSequence = seqTop;
            mSideSequence = seqSide;
            mMatchScore = matchScore;
            mMismatchScore = mismatchScore;
            mGapScore = gapScore;

            var width = mTopSequence.length + 1;
            var height = mSideSequence.length + 1;

            mPathTable = new Array(width).fill(0).map(() => new Array(height).fill(0));
            mCellMap = {};

            for (var i = 0; i < width; i++) {
                for (var j = 0; j < height; j++) {
                    if (i === 0 || j === 0) {
                        mPathTable[i][j] = 0;
                        mCellMap[i + "_" + j] = {
                            'winningScore': mPathTable[i][j]
                        };
                        continue;
                    }

                    var isMatch = mTopSequence[i - 1] === mSideSequence[j - 1];
                    var comparisonScore = isMatch ? matchScore : mismatchScore;

                    var moveUpScore = mPathTable[i][j - 1] + gapScore;
                    var moveSdScore = mPathTable[i - 1][j] + gapScore;
                    var moveDgScore = mPathTable[i - 1][j - 1] + comparisonScore;
                    mPathTable[i][j] = Math.max(0, moveUpScore, moveSdScore, moveDgScore);

                    var direction = [];
                    if (mPathTable[i][j] === moveDgScore) {
                        direction.push('d');
                    }

                    if (mPathTable[i][j] === moveUpScore) {
                        direction.push('u');
                    }

                    if (mPathTable[i][j] === moveSdScore) {
                        direction.push('s');
                    }

                    mCellMap[i + "_" + j] = {
                        'sideScoreText': mPathTable[i - 1][j] + " + " + gapScore + " (The Gap score) = " + moveSdScore,
                        'upScoreText': mPathTable[i][j - 1] + " + " + gapScore + " (The Gap score) = " + moveUpScore,
                        'diagonalScoreText': mPathTable[i - 1][j - 1] + " + " + comparisonScore + " (Due to a " + (isMatch ? "match" : "mismatch") + " between " + mTopSequence[i - 1] + " & " + mSideSequence[j - 1] + ") = " + moveDgScore,
                        'sideScore': moveSdScore,
                        'upScore': moveUpScore,
                        'diagonalScore': moveDgScore,
                        'winningScore': mPathTable[i][j],
                        'direction': direction
                    };
                }
            }

            constructGrid();
        }
    };

    return mSelf;
}());
